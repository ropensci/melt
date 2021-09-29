#include "utils_gbd.h"

// [[Rcpp::export]]
Rcpp::List test(const Eigen::MatrixXd& x,
                const Eigen::MatrixXd& c,
                const Eigen::MatrixXd& lhs,
                const Eigen::VectorXd& rhs,
                const double threshold,
                const int maxit = 1e4,
                const double abstol = 1e-8) {
  // check lhs & rhs
  if (lhs.rows() > lhs.cols()) {
    Rcpp::stop("nrow(lhs) must not exceed ncol(lhs).");
  }
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(lhs);
  if (lu_decomp.rank() != lhs.rows()) {
    Rcpp::stop("lhs must have full rank.");
  }
  if (lhs.rows() != rhs.rows()) {
    Rcpp::stop("dimensions of lhs and rhs do not match.");
  }

  // global minimizer
  const Eigen::VectorXd theta_hat =
    x.array().colwise().sum() / c.array().colwise().sum();

  minEL el =
    test_gbd_EL(theta_hat, x, c, lhs, rhs, threshold, maxit, abstol);

  Rcpp::List result;
  result["theta"] = el.theta;
  result["lambda"] = el.lambda;
  result["n2logLR"] = 2 * el.nlogLR;
  result["convergence"] = el.convergence;
  result["iterations"] = el.iterations;
  result.attr("class") = Rcpp::CharacterVector({"test", "melt"});
  return result;

}
