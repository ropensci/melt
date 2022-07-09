#include "test_multiple_hypotheses_utils.h"
#include "utils.h"
#include "EL.h"
#include <RcppEigen.h>
#include <string>
#include <vector>

// [[Rcpp::export]]
Rcpp::List test_multiple_hypotheses(const double alpha,
                                    const Eigen::Map<Eigen::VectorXi> &q,
                                    const int m,
                                    const int M,
                                    const std::string method,
                                    const Eigen::Map<Eigen::VectorXd> &est,
                                    const Eigen::Map<Eigen::MatrixXd> &x,
                                    const Eigen::Map<Eigen::VectorXd> &rhs,
                                    const Eigen::Map<Eigen::MatrixXd> &lhs,
                                    const int maxit,
                                    const int maxit_l,
                                    const double tol,
                                    const double tol_l,
                                    const Rcpp::Nullable<double> step,
                                    const Rcpp::Nullable<double> th,
                                    const Eigen::Map<Eigen::ArrayXd> &w)
{
  const double gamma = set_step(x.rows(), step);
  // 1. compute test statistics
  std::vector<double> test_statistic(m);
  for (int j = 0; j < m; ++j)
  {
    const Eigen::VectorXd r = rhs.middleRows(q(j), q(j + 1) - q(j));
    const Eigen::MatrixXd l = lhs.middleRows(q(j), q(j + 1) - q(j));
    const double test_th = set_threshold(l.rows(), th);
    const CEL el(method, est, x, l, r, maxit, maxit_l, tol, tol_l, gamma,
                 test_th, w);
    test_statistic[j] = 2.0 * el.nllr;
  }
  // 2. compute critical value
  const Eigen::MatrixXd s = shat(method, est, x);
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
  const Eigen::MatrixXd sqrt_s = es.operatorSqrt();
  // ahat matrices
  const int p = est.size();
  const Eigen::MatrixXd h = dg0_inv(method, x, est);
  Eigen::MatrixXd amat(p, p * m);
  for (int j = 0; j < m; ++j)
  {
    const Eigen::MatrixXd l = lhs.middleRows(q(j), q(j + 1) - q(j));
    amat.middleCols(j * p, p) = ahat(l, h, s);
  }
  // Monte Carlo simulation
  std::vector<double> max_statistic(M);
  // #pragma omp parallel for num_threads(nthreads)
  for (int b = 0; b < M; ++b)
  {
    const Eigen::RowVectorXd u = rmvn(sqrt_s);
    Eigen::MatrixXd tmp = u * amat;
    tmp.resize(p, m);
    max_statistic[b] = (u * tmp).maxCoeff();
  }
  const double cv = compute_quantile(Rcpp::wrap(max_statistic), 1 - alpha);
  // 3. compute adjusted p-values
  std::vector<double> adj_pval(m);
  for (int j = 0; j < m; ++j)
  {
    adj_pval[j] =
        count_if(
            max_statistic.begin(), max_statistic.end(),
            [test_statistic, j](double x)
            { return (x > test_statistic[j]); }) /
        static_cast<double>(M);
  }
  Rcpp::List result = Rcpp::List::create(
      Rcpp::Named("statistic") = Rcpp::wrap(test_statistic),
      Rcpp::Named("cv") = cv,
      Rcpp::Named("pval") = adj_pval);
  return result;
}
