#include "testMultipleHypotheses_utils.h"
#include "utils.h"
#include "EL.h"
#include <RcppEigen.h>
#include <string>
#include <vector>

// [[Rcpp::export]]
Rcpp::List testMultipleHypotheses(const Eigen::Map<Eigen::VectorXi>& q,
                                  const int m,
                                  const int B,
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
  // const Rcpp::NumericVector out = Rcpp::wrap(test_statistic);
  // return out;



  // sample covariance with plug-in estimate
  const Eigen::MatrixXd s = cov(method, est, x);
  // sqrt of sample covariance
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
  const Eigen::MatrixXd sqrt_s = es.operatorSqrt();

  // ahat matrices
  const int p = est.size();
  const Eigen::MatrixXd h = dg0_inv(method, x);
  Eigen::MatrixXd amat(p, p * m);
  for (int j = 0; j < m; ++j) {
    const Eigen::MatrixXd l = lhs.middleRows(q(j), q(j + 1) - q(j));
    amat.middleCols(j * p, p) = ahat(l, h, s);
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
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("tmp") = Rcpp::wrap(test_statistic),
    Rcpp::Named("tmp2") = Rcpp::wrap(max_statistic));
  return result;
}
