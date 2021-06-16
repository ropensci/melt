#include "utils.h"

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
Rcpp::List el_mean(arma::rowvec theta, arma::mat x,
               int maxit = 100, double abstol = 1e-8) {
  if (arma::rank(x) != x.n_cols) {
    Rcpp::stop("Design matrix x must have full rank.");
  }
  // ncol must be same as p
  const int n = x.n_rows;
  const int p = theta.n_elem;

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
    f0 = -arma::sum(plog(arg));
    // J matrix & y vector
    arma::vec v1 = arma::sqrt(-d2plog(arg));
    arma::vec v2 = dplog(arg);
    J = g.each_col() % v1;
    y = v2 / v1;
    // update lambda by NR method with least square & step halving
    arma::qr_econ(Q, R, J);
    lc = l + arma::solve(R, trans(Q) * y);
    double alpha = 1;
    while (-arma::sum(plog(1 + g * lc)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * solve(R, trans(Q) * y);
    }
    // update function value
    l = lc;
    arg = 1 + g * l;
    f1 = -arma::sum(plog(arg));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      arma::vec v1 = arma::sqrt(-d2plog(arg));
      arma::vec v2 = dplog(arg);
      J = g.each_col() % v1;
      y = v2 / v1;
      convergence = true;
    } else {
      ++iterations;
      if (iterations == maxit) {
        break;
      }
    }
  }

  // result
  return Rcpp::List::create(
    Rcpp::Named("nlogLR") = -f1,
    Rcpp::Named("lambda") = Rcpp::NumericVector(l.begin(), l.end()),
    Rcpp::Named("grad") = Rcpp::as<std::vector<double>>(Rcpp::wrap(-trans(J) * y)),
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("convergence") = convergence);
}

//' Empirical likelihood test for mean
//'
//' Compute empirical likelihood for mean
//'
//' @param theta a vector of parameters to be tested.
//' @param x a matrix or vector of data. Each row is an observation vector.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 50.
//' @export
// [[Rcpp::export]]
Rcpp::List el_mean2(const arma::vec& theta,
                    const arma::mat& x,
                    const int maxit = 100,
                    const double abstol = 1e-8) {
  if (arma::rank(x) != x.n_cols) {
    Rcpp::stop("Design matrix x must have full rank.");
  }
  // estimating function for mean parameters
  // const arma::mat g = g_mean(theta, x);
  // compute EL
  const EL result = getEL(g_mean(theta, x), maxit, abstol);

  return Rcpp::List::create(
    Rcpp::Named("nlogLR") = result.nlogLR,
    Rcpp::Named("lambda") =
      Rcpp::NumericVector(result.lambda.begin(), result.lambda.end()),
    Rcpp::Named("gradient") =
      Rcpp::NumericVector(result.gradient.begin(), result.gradient.end()),
    Rcpp::Named("iterations") = result.iterations,
    Rcpp::Named("convergence") = result.convergence);
}



