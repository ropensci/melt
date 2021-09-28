#include "EL.h"

//' Empirical likelihood test for mean
//'
//' Compute empirical likelihood for mean
//'
//' @param theta Numeric vector of parameters to be tested.
//' @param x Numeric matrix or vector of data. If \code{x} is a matrix, each row corresponds to an observation.
//' @param maxit Maximum number of iterations for optimization. Defaults to 50.
//' @param abstol Absolute convergence tolerance for optimization. Defaults to 1e-08.
//'
//' @return A list with class \code{elmulttest}.
//'
//' @examples
//' ## scalar mean
//' theta <- 0
//' x <- rnorm(100)
//' el_mean(theta, x)
//'
//' ## vector mean
//' x <- matrix(rnorm(100), ncol = 2)
//' theta <- c(0, 0)
//' el_mean(theta, x)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List el_mean(const Eigen::Map<Eigen::VectorXd>& theta,
                   const Eigen::Map<Eigen::MatrixXd>& x,
                   const int maxit = 50,
                   const double abstol = 1e-8) {
  if (theta.size() != x.cols()) {
    Rcpp::stop("dimensions of theta ans x do not match.");
  }
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != x.cols()) {
    Rcpp::stop("design matrix x must have full rank.");
  }
  // compute EL
  // const EL result = getEL(x.rowwise() - theta.transpose(), maxit, abstol);
  const EL2 result(x.rowwise() - theta.transpose(), x.cols() * 100, maxit, abstol);

  return Rcpp::List::create(
    Rcpp::Named("n2logLR") = 2 * result.nlogLR,
    Rcpp::Named("lambda") = result.lambda,
    Rcpp::Named("iterations") = result.iterations,
    Rcpp::Named("convergence") = result.convergence);
}

