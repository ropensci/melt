#include "EL.h"

// [[Rcpp::export]]
Eigen::MatrixXd confint_(
    const std::string method,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const double cutoff,
    const Rcpp::IntegerVector& idx,
    const int maxit,
    const int maxit_l,
    const double tol,
    const double tol_l,
    const Rcpp::Nullable<double> step,
    const Rcpp::Nullable<double> th,
    const int nthreads,
    const Eigen::Map<Eigen::ArrayXd>& w)
{
  // parameter dimension
  const int p = par0.size();
  // number of confidence intervals
  const int n = idx.size();
  // confidence intervals (vector)
  std::vector<double> ci_vec;
  ci_vec.reserve(2 * n);
  // step size
  const double gamma = step_nloglr(x.rows(), step);
  // test threshold
  const double test_th = th_nloglr(1, th);
  // #pragma omp parallel for num_threads(nthreads)
  for (int j : idx) {
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(j - 1) = 1.0;

    // lower endpoint
    double lower_ub = par0[j - 1];
    double lower_lb = par0[j - 1] - 1.0 / std::log(x.rows());
    // lower bound for lower endpoint
    while (2.0 * MINEL(method, par0, x, lhs,
                       Eigen::Matrix<double, 1, 1>(lower_lb), maxit, maxit_l,
                       tol, tol_l, gamma, test_th, w).nllr <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= 1.0 / std::log(x.rows());
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > tol) {
      if (2.0 * MINEL(method, par0, x, lhs,
                      Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2.0),
                      maxit, maxit_l, tol, tol_l, gamma, test_th, w).nllr >
            cutoff) {
        lower_lb = (lower_lb + lower_ub) / 2.0;
      } else {
        lower_ub = (lower_lb + lower_ub) / 2.0;
      }
    }
    ci_vec.push_back(lower_ub);

    // upper endpoint
    double upper_lb = par0[j - 1];
    double upper_ub = par0[j - 1] + 1.0 / std::log(x.rows());
    // upper bound for upper endpoint
    while (2.0 * MINEL(method, par0, x, lhs,
                       Eigen::Matrix<double, 1, 1>(upper_ub), maxit, maxit_l,
                       tol, tol_l, gamma, test_th, w).nllr <= cutoff) {
      upper_lb = upper_ub;
      upper_ub += 1.0 / std::log(x.rows());
    }
    // approximate upper bound by numerical search
    while (upper_ub - upper_lb > tol) {
      if (2.0 * MINEL(method, par0, x, lhs,
                      Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2.0),
                      maxit, maxit_l, tol, tol_l, gamma, test_th, w).nllr >
            cutoff) {
        upper_ub = (upper_lb + upper_ub) / 2.0;
      } else {
        upper_lb = (upper_lb + upper_ub) / 2.0;
      }
    }
    ci_vec.push_back(upper_lb);
  }

  // transform into a matrix
  Eigen::MatrixXd ci = Eigen::Map<Eigen::MatrixXd>(ci_vec.data(), 2, n);
  return ci.transpose();
}
