#include "utils.h"
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
Rcpp::List el_mean(const Eigen::Map<Eigen::VectorXd>& theta,
                   const Eigen::Map<Eigen::MatrixXd>& x,
                   const int maxit = 100,
                   const double abstol = 1e-8) {
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  // // auto rank = Eigen::ColPivHouseholderQR< Eigen::MatrixXd >::rank(x);
  if (lu_decomp.rank() != x.cols()) {
    Rcpp::stop("Design matrix x must have full rank.");
  }
  // compute EL
  const EL& result = getEL(x.rowwise() - theta.transpose(), maxit, abstol);

  // const TEST result(x.rowwise() - theta.transpose(), maxit, abstol);

  return Rcpp::List::create(
    Rcpp::Named("nlogLR") = result.nlogLR,
    Rcpp::Named("lambda") = result.lambda,
    Rcpp::Named("iterations") = result.iterations,
    Rcpp::Named("convergence") = result.convergence);
}

