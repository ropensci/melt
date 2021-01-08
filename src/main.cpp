#include <RcppArmadillo.h>
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


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
  arma::mat g = x;
  g.each_row() -= theta;

  // minimization
  arma::vec l; l.zeros(p);
  arma::vec lc;
  arma::vec arg = 1 + g * l;
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
    f0 = -sum(plogcpp(wrap(arg), 1 / n));
    // J matrix & y vector
    arma::vec v1(Rcpp::sqrt(-d2plogcpp(wrap(arg), 1 / n)));
    arma::vec v2(dplogcpp(wrap(arg), 1 / n));
    J = g.each_col() % v1;
    y = v2 / v1;
    // update lambda by NR method with least square & step halving
    arma::qr_econ(Q, R, J);
    lc = l + arma::solve(R, trans(Q) * y);
    double alpha = 1;
    while(-sum(plogcpp(wrap(1 + g * lc), 1 / n)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * solve(R, trans(Q) * y);
    }
    // update function value
    l = lc;
    arg = 1 + g * l;
    f1 = -sum(plogcpp(wrap(arg), 1 / n));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      arma::vec v1(Rcpp::sqrt(-d2plogcpp(wrap(arg), 1 / n)));
      arma::vec v2(dplogcpp(wrap(arg), 1 / n));
      J = g.each_col() % v1;
      y = v2 / v1;
      convergence = true;
    } else {
      iterations++;
      if(iterations == maxit) {
        break;
      }
    }
  }

  // result
  return Rcpp::List::create(
    Rcpp::Named("nlogLR") = -f1,
    Rcpp::Named("lambda") = NumericVector(l.begin(), l.end()),
    Rcpp::Named("grad") = as<std::vector<double>>(wrap(-trans(J) * y)),
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("convergence") = convergence);
}


//' Two sample test for equal mean
//'
//' Two sample test for equal mean
//'
//' @param x a vector of data for one  group.
//' @param y a vector of data for the other  group.
//' @param b a momentum parameter for minimization. Defaults to .1.
//' @param alpha an optional step size. Defaults to 1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @examples
//' x <- rnorm(100)
//' y <- rnorm(100)
//' test2sample2_cpp(x, y)
//'
//' @export
// [[Rcpp::export]]
List test2sample2_cpp(NumericVector x, NumericVector y, double b = .9, double alpha = 1,
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
  unsigned int nx = x.size();
  unsigned int ny = y.size();
  unsigned int N = std::max(nx, ny);
  unsigned int iterations{};
  int convergence{};
  double v{};
  double lx;
  double ly;
  Rcpp::NumericVector grad;
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
  return Rcpp::List::create(
    Rcpp::Named("par") = par,
    Rcpp::Named("nlogLR") = sum(Rcpp::log(1 + (x - par) * lx)) +
      sum(Rcpp::log(1 + (y - par) * ly)),
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("convergence") = convergence
  );
}


//' Two sample test for equal mean
//'
//' Two sample test for equal mean
//'
//' @param x a vector of data for one  group.
//' @param y a vector of data for the other  group.
//' @param b a momentum parameter for minimization. Defaults to .1.
//' @param alpha an optional step size. Defaults to 1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @examples
//' x <- rnorm(100)
//' y <- rnorm(100)
//' test2sample777_cpp(x, y)
//'
//' @export
// [[Rcpp::export]]
List test2sample777_cpp(NumericVector x, NumericVector y, double b = .9, double alpha = 1,
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
  unsigned int nx = x.size();
  unsigned int ny = y.size();
  unsigned int N = std::max(nx, ny);
  unsigned int iterations{};
  int convergence{};
  double v{};
  double lx;
  double ly;
  Rcpp::NumericVector grad;
  double d;

  // linearization
  double par0 = par;
  double lx0 = elMeancpp({par}, x)["lambda"];
  double ly0 = elMeancpp({par}, y)["lambda"];
  double mx = -sum(1 / pow(1 + lx0 * (x - par), 2)) /
    sum(pow(x - par, 2) / pow(1 + lx0 * (x - par), 2));
  double my = -sum(1 / pow(1 + ly0 * (y - par) , 2)) /
    sum(pow(y - par ,2) / pow(1 + ly0 * (y - par) , 2));

  // minimization
  while (convergence == 0) {
    // lambda update
    lx = mx * (par - par0) + lx0;
    ly = my * (par - par0) + ly0;
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
  return Rcpp::List::create(
    Rcpp::Named("par") = par,
    Rcpp::Named("nlogLR") = sum(Rcpp::log(1 + (x - par) * lx)) +
      sum(Rcpp::log(1 + (y - par) * ly)),
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("convergence") = convergence
  );
}
