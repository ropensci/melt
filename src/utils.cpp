#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec plogcpp(arma::vec x, double threshold) {
  arma::vec result(x.n_rows);
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
// [[Rcpp::export]]
arma::vec dplogcpp(arma::vec x, double threshold) {
  arma::vec result(x.n_rows);
  for (uword i = 0; i < result.n_rows; i++) {
    if(x(i) >= threshold) {
      result(i) = pow(x(i), -1);
    } else {
      result(i) = 2 * pow(threshold, -1) - x(i) * pow(threshold, -2);
    }
  }
  return result;
}
// [[Rcpp::export]]
arma::vec d2plogcpp(arma::vec x, double threshold) {
  arma::vec result(x.n_rows);
  for (uword i = 0; i < result.n_rows; i++) {
    if(x(i) >= threshold) {
      result(i) = -pow(x(i), -2);
    } else {
      result(i) = -pow(threshold, -2);
    }
  }
  return result;
}
// [[Rcpp::export]]
double dfp1dcpp(arma::vec gr) {
  arma::mat LHS = {{1, 0, 1}, {0, 1, -1}, {1, -1, 0}};
  arma::vec RHS = zeros<vec>(3);
  RHS.subvec(0, 1) = -gr;
  arma::vec sol = solve(LHS, RHS);
  return sol(0);
}


// [[Rcpp::export]]
List elMeancpp(arma::rowvec theta, arma::mat x,
               int maxit = 100, double abstol = 1e-8) {
  // ncol must be same as p
  int n = x.n_rows;
  int p = theta.n_elem;
  // estimating equation
  arma::mat g = x;
  g.each_row() -= theta;

  // minimization
  arma::vec l; l.zeros(p);
  arma::vec lc;
  arma::vec y;
  arma::mat J;
  arma::mat Q;
  arma::mat R;
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


// [[Rcpp::export]]
List test2sample_cpp(arma::vec x, arma::vec y, double b = .1,
                     int maxit = 1000, double abstol = 1e-8) {
  List result;
  double ub = std::min(x.max(), y.max());
  // colvec ub2 = conv_to< colvec >::from(ub);
  double lb = std::max(x.min(), y.min());
  if(ub <= lb) {
    result["nlogLR"] = std::numeric_limits<double>::infinity();
    result["convergence"] = -1;
    return result;
  }

  // initialization
  double par = (lb + ub) / 2;
  int nx = x.n_elem;
  int ny = y.n_elem;
  int N = nx + ny;
  double alpha = pow(N, -1);
  int iterations = 0;
  int convergence = 0;
  double v = 0;
  double lx;
  double ly;
  arma::vec grad;
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
      v = -v;
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
        v = v / 2;
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
