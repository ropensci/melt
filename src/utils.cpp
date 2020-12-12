#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

vec plogcpp(vec x, double threshold) {
  vec result(x.n_rows);
  for (uword i = 0; i < result.n_rows; i++) {
    if(x(i) >= threshold) {
      result(i) = log(x(i));
    } else {
      result(i) = log(threshold) - 1.5 + 2 * pow(threshold, -1) * x(i) -
        pow(threshold, -2) * pow(x(i), 2) / 2;
    }
  }
  return result;
}
vec dplogcpp(vec x, double threshold) {
  vec result(x.n_rows);
  for (uword i = 0; i < result.n_rows; i++) {
    if(x(i) >= threshold) {
      result(i) = pow(x(i), -1);
    } else {
      result(i) = 2 * pow(threshold, -1) - x(i) * pow(threshold, -2);
    }
  }
  return result;
}
vec d2plogcpp(vec x, double threshold) {
  vec result(x.n_rows);
  for (uword i = 0; i < result.n_rows; i++) {
    if(x(i) >= threshold) {
      result(i) = -pow(x(i), -2);
    } else {
      result(i) = -pow(threshold, -2);
    }
  }
  return result;
}
double dfp1dcpp(vec gr) {
  mat LHS = {{1, 0, 1}, {0, 1, -1}, {1, -1, 0}};
  vec RHS = zeros<vec>(3);
  RHS.subvec(0, 1) = -gr;
  vec sol = solve(LHS, RHS);
  return sol(0);
}



//' Compute empirical likelihood for mean
//'
//' Compute empirical likelihood for mean
//'
//' @param theta a vector of parameters to be tested.
//' @param x a matrix or vector of data. Each row is an observation vector.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 50.
//' @export
// [[Rcpp::export]]
List elMeancpp(arma::rowvec theta, arma::mat x,
               int maxit = 100, double abstol = 1e-8) {
  // ncol must be same as p
  int n = x.n_rows;
  int p = theta.n_elem;
  // estimating equation
  mat g = x;
  g.each_row() -= theta;

  // minimization
  vec l; l.zeros(p);
  vec lc;
  vec y;
  mat J;
  mat Q;
  mat R;
  double f0;
  double f1;
  int iterations = 0;
  bool convergence = false;
  while (convergence == false) {
    // function evaluation(initial)
    f0 = -sum(plogcpp(1 + g * l, 1 / n));
    // J matrix for least square
    J = g.each_col() % sqrt(-d2plogcpp(1 + g * l, 1 / n));
    // Y vector for least square
    y = dplogcpp(1 + g * l, 1 / n) / sqrt(-d2plogcpp(1 + g * l, 1 / n));
    // update lambda by NR method with least square & step halving
    qr_econ(Q, R, J);
    lc = l + solve(R, trans(Q) * y);
    double alpha = 1;
    while (-sum(plogcpp(1 + g * lc, 1 / n)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * solve(R, trans(Q) * y);
    }
    // update function value
    l = lc;
    f1 = -sum(plogcpp(1 + g * l, 1 / n));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      convergence = true;
    } else {
      iterations++;
      if(iterations == maxit) {
        break;
      }
    }
  }

  // result
  List result;
  result["nlogLR"] = sum(plogcpp(1 + g * l, 1 / n));
  result["lambda"] = NumericVector(l.begin(), l.end());
  result["grad"] = sum(g.each_col() % dplogcpp(1 + g *l, 1 / n));
  result["iterations"] = iterations;
  result["convergence"] = convergence;
  return result;
}


//' Two sample test for equal mean
//'
//' Two sample test for equal mean
//'
//' @param x a vector of data for one  group.
//' @param y a vector of data for the other  group.
//' @param b a momentum parameter for minimization. Defaults to .1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @examples
//' x <- rnorm(100)
//' y <- rnorm(100)
//' test2sample2(x, y)
//'
//' @export
// [[Rcpp::export]]
List test2sample_cpp(arma::vec x, arma::vec y, double b = .9,
                     unsigned int maxit = 1000, double abstol = 1e-8) {
  List result;
  double ub = std::min(x.max(), y.max());
  double lb = std::max(x.min(), y.min());
  if(ub <= lb) {
    result["nlogLR"] = datum::inf;
    result["convergence"] = -1;
    return result;
  }

  // initialization
  double par = (lb + ub) / 2;
  unsigned int nx = x.n_elem;
  unsigned int ny = y.n_elem;
  double alpha = (ub - lb) / (nx + ny);
  unsigned int iterations{};
  int convergence{};
  double v{};
  double lx;
  double ly;
  vec grad;
  double d;

  // minimization
  while (convergence == 0) {
    // lambda update
    lx = elMeancpp({par}, x)["lambda"];
    ly = elMeancpp({par}, y)["lambda"];
    // gradient
    grad = {sum(dplogcpp(1 + lx * (x - par), 1 / nx)) * (-lx),
      sum(dplogcpp(1 + ly * (y - par), 1 / ny)) * (-ly)};
    // direction
    d = dfp1dcpp(grad);
    // direction change reverts momentum
    if (sign(d) != sign(v)) {
      v = 0; alpha = alpha / 2;
    }
    // lb, ub update
    if (sign(d) > 0) {
      lb = par;
    } else {
      ub = par;
    }
    // convergence check & parameter update
    if (std::abs(d * sum(grad)) < abstol || ub - lb < abstol) {
      convergence = 1;
    } else {
      iterations++;
      // step halving to satisfy convex hull constraint
      v = b * v + d;
      while (par + alpha * v <= lb || par + alpha * v >= ub) {
        alpha = alpha / 2;
      }
      par = par + alpha * v;
      if(iterations == maxit) {
        break;
      }
    }
  }

  // result
  result["par"] = par;
  double xnlogLR = elMeancpp({par}, x)["nlogLR"];
  double ynlogLR = elMeancpp({par}, y)["nlogLR"];
  result["nlogLR"] = xnlogLR  + ynlogLR;
  result["iterations"] = iterations;
  result["convergence"] = convergence;
  return result;
}


//' Two sample test for equal mean
//'
//' Two sample test for equal mean
//'
//' @param x a vector of data for one  group.
//' @param y a vector of data for the other  group.
//' @param b a momentum parameter for minimization. Defaults to .1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @examples
//' x <- rnorm(100)
//' y <- rnorm(100)
//' test2sample2(x, y)
//'
//' @export
// [[Rcpp::export]]
List test2sample2_cpp(arma::vec x, arma::vec y, double b = .9, double alpha = 1,
                      unsigned int maxit = 1000, double abstol = 1e-8) {
  List result;
  vec sample_mean = {mean(x), mean(y)};
  vec ub_vec = {max(x), max(y), max(sample_mean)};
  vec lb_vec = {min(x), min(y), min(sample_mean)};
  double ub = min(ub_vec);
  double lb = max(lb_vec);
  if(ub <= lb) {
    result["nlogLR"] = datum::inf;;
    result["convergence"] = -1;
    return result;
  }

  // initialization
  double par = (lb + ub) / 2;
  unsigned int nx = x.n_elem;
  unsigned int ny = y.n_elem;
  unsigned int N = std::max(nx, ny);
  unsigned int iterations{};
  int convergence{};
  double v{};
  double lx;
  double ly;
  vec grad;
  double d;

  // minimization
  while (convergence == 0) {
    // lambda update
    lx = elMeancpp({par}, x)["lambda"];
    ly = elMeancpp({par}, y)["lambda"];
    // gradient
    grad = {sum(dplogcpp(1 + lx * (x - par), 1 / nx)) * (-lx) / N,
            sum(dplogcpp(1 + ly * (y - par), 1 / ny)) * (-ly) / N};
    // direction
    d = dfp1dcpp(grad);
    // direction change reverts momentum
    if (sign(d) != sign(v)) {
      v = 0; alpha = alpha / 2;
    }
    // lb, ub update
    if (sign(d) > 0) {
      lb = par;
    } else {
      ub = par;
    }
    // convergence check & parameter update
    if ((std::abs(d * sum(grad)) < abstol || ub - lb < abstol) && iterations > 0) {
      convergence = 1;
    } else {
      iterations++;
      // step halving to satisfy convex hull constraint
      v = b * v + d;
      while (par + alpha * v <= lb || par + alpha * v >= ub) {
        alpha = alpha / 2;
      }
      par = par + alpha * v;
      if(iterations == maxit) {
        break;
      }
    }
  }

  // result
  result["par"] = par;
  result["nlogLR"] = sum(log(1 + (x - par) * lx)) +
    sum(log(1 + (y - par) * ly));
  result["iterations"] = iterations;
  result["convergence"] = convergence;
  return result;
}
