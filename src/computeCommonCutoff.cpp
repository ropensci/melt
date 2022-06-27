#include "elmt_utils.h"
#include <RcppEigen.h>
#include <string>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector computeCommonCutoff(const std::string method,
                                        const Eigen::Map<Eigen::VectorXd>& est,
                                        const Eigen::Map<Eigen::MatrixXd>& x,
                                        const Eigen::Map<Eigen::MatrixXd>& lhs,
                                        const Eigen::Map<Eigen::VectorXi>& q,
                                        const int m,
                                        const int B) {
  // sample covariance with plug-in estimate
  const Eigen::MatrixXd s = cov(method, est, x);
  // sqrt of sample covariance
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
  const Eigen::MatrixXd sqrt_s = es.operatorSqrt();

  // ahat matrices
  const int p = est.size();
  const Eigen::MatrixXd w = dg0_inv(method, x);
  Eigen::MatrixXd amat(p, p * m);
  for (int j = 0; j < m; ++j) {
    const Eigen::MatrixXd l = lhs.middleRows(q(j), q(j + 1) - q(j));
    amat.middleCols(j * p, p) = ahat(l, w, s);
  }

  // 4. random monte carlo
  std::vector<double> max_statistic(B);
  // #pragma omp parallel for num_threads(nthreads)
  for (int b = 0; b < B; ++b) {
    const Eigen::RowVectorXd u = rmvn(sqrt_s);
    Eigen::MatrixXd tmp = u * amat;
    tmp.resize(p, m);
    max_statistic[b] = (u * tmp).maxCoeff();
  }

  // critical value can be computed independent of the statistics
  const Rcpp::NumericVector out = Rcpp::wrap(max_statistic);
  return out;
}

// // [[Rcpp::export]]
// Rcpp::List elmt_(const std::string method,
//                  const Eigen::Map<Eigen::VectorXd>& est,
//                  const Eigen::Map<Eigen::MatrixXd>& x,
//                  const Eigen::Map<Eigen::VectorXd>& rhs,
//                  const Eigen::Map<Eigen::MatrixXd>& lhs,
//                  const Eigen::Map<Eigen::VectorXi>& q,
//                  const int m,
//                  const int B,
//                  const int maxit,
//                  const int maxit_l,
//                  const double tol,
//                  const double tol_l,
//                  const Rcpp::Nullable<double> step,
//                  const Rcpp::Nullable<double> th,
//                  const Eigen::Map<Eigen::ArrayXd>& wt) {
//   const double gamma = step_nloglr(x.rows(), step);
//   // test statistics
//   // 1. compute the vector of test statistics
//   std::vector<double> test_statistic(m);
//   for (int j = 0; j < m; ++j) {
//     const Eigen::MatrixXd l = lhs.middleRows(q(j), q(j + 1) - q(j));
//     const Eigen::VectorXd r = rhs.middleRows(q(j), q(j + 1) - q(j));
//     const double test_th = th_nloglr(l.rows(), th);
//     const MINEL el(method, est, x, l, r, maxit, maxit_l, tol, tol_l, gamma,
//                    test_th, wt);
//     test_statistic[j] = 2.0 * el.nllr;
//   }
//
//   // 2. compute sample covariance with estimate(given in R)
//   const Eigen::MatrixXd w = dg0_inv(method, x);
//   // sample covariance
//   const Eigen::MatrixXd s = cov(method, est, x);
//   // sqrt of sample covariance
//   const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
//   const Eigen::MatrixXd sqrt_s = es.operatorSqrt();
//
//   // 3. prepare ahat matrices
//   const int p = est.size();
//   Eigen::MatrixXd amat(p, p * m);
//   for (int j = 0; j < m; ++j) {
//     const Eigen::MatrixXd l = lhs.middleRows(q(j), q(j + 1) - q(j));
//     amat.middleCols(j * p, p) = ahat(l, w, s);
//   }
//
//   // 4. random monte carlo
//   std::vector<double> max_statistic(B);
//   // #pragma omp parallel for num_threads(nthreads)
//   for (int b = 0; b < B; ++b) {
//     const Eigen::RowVectorXd u = rmvn(sqrt_s);
//     Eigen::MatrixXd tmp = u * amat;
//     tmp.resize(p, m);
//     max_statistic[b] = (u * tmp).maxCoeff();
//   }
//
//   // critical value can be computed independent of the statistics
//   Rcpp::List result = Rcpp::List::create(
//     Rcpp::Named("tmp") = Rcpp::wrap(test_statistic),
//     Rcpp::Named("tmp2") = Rcpp::wrap(max_statistic));
//   return result;
// }
