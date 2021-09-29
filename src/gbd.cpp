#include "utils_gbd.h"

// [[Rcpp::export]]
Rcpp::List test_gbd(const Eigen::MatrixXd& x,
                    const Eigen::MatrixXd& c,
                    const Eigen::MatrixXd& lhs,
                    const Eigen::VectorXd& rhs,
                    const bool approx = false,
                    const int maxit = 1000,
                    const double abstol = 1e-8) {
  /// initialization ///
  // if (arma::rank(L) != lhs.rows()) {
  //   Rcpp::stop("Hypothesis matrix lhs must have full rank.");
  // }
  if (lhs.rows() != rhs.rows()) {
    Rcpp::stop("Dimensions of L and rhs do not match.");
  }
  minEL result =
    test_gbd_EL(x.array().colwise().sum() / c.array().colwise().sum(),
                x, c, lhs, rhs, maxit, abstol);

  return Rcpp::List::create(
    Rcpp::Named("theta") = result.theta,
    Rcpp::Named("lambda") = result.lambda,
    Rcpp::Named("nlogLR") = result.nlogLR,
    Rcpp::Named("iterations") = result.iterations,
    Rcpp::Named("convergence") = result.convergence);
}
