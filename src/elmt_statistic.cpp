#include "utils.h"
#include "EL.h"
#include <RcppEigen.h>
#include <string>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector elmt_statistic_(const Eigen::Map<Eigen::VectorXi>& q,
                                    const int m,
                                    const std::string method,
                                    const Eigen::Map<Eigen::VectorXd>& est,
                                    const Eigen::Map<Eigen::MatrixXd>& x,
                                    const Eigen::Map<Eigen::VectorXd>& rhs,
                                    const Eigen::Map<Eigen::MatrixXd>& lhs,
                                    const int maxit,
                                    const int maxit_l,
                                    const double tol,
                                    const double tol_l,
                                    const Rcpp::Nullable<double> step,
                                    const Rcpp::Nullable<double> th,
                                    const Eigen::Map<Eigen::ArrayXd>& w) {
  const double gamma = step_nloglr(x.rows(), step);
  std::vector<double> test_statistic(m);
  for (int j = 0; j < m; ++j) {
    const Eigen::VectorXd r = rhs.middleRows(q(j), q(j + 1) - q(j));
    const Eigen::MatrixXd l = lhs.middleRows(q(j), q(j + 1) - q(j));
    const double test_th = th_nloglr(l.rows(), th);
    const MINEL el(method, est, x, l, r, maxit, maxit_l, tol, tol_l, gamma,
                   test_th, w);
    test_statistic[j] = 2.0 * el.nllr;
  }
  const Rcpp::NumericVector out = Rcpp::wrap(test_statistic);
  return out;
}
