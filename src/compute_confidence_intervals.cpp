#include "utils.h"
#include "EL.h"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string>

// [[Rcpp::export]]
Eigen::MatrixXd compute_confidence_intervals(
    const std::string method,
    const Eigen::Map<Eigen::VectorXd> &par0,
    const Eigen::Map<Eigen::MatrixXd> &x,
    const double cutoff,
    const Rcpp::IntegerVector &idx,
    const int maxit,
    const int maxit_l,
    const double tol,
    const double tol_l,
    const Rcpp::Nullable<double> step,
    const Rcpp::Nullable<double> th,
    const int nthreads,
    const Eigen::Map<Eigen::ArrayXd> &w)
{
  // parameter length
  const int p = par0.size();
  // number of confidence intervals
  const int n = idx.size();
  // matrix of confidence intervals
  Eigen::MatrixXd ci(n, 2);

  // step size
  const double gamma = set_step(x.rows(), step);
  // test threshold
  const double test_th = set_threshold(1, th);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nthreads)
  #endif
  for (int i = 0; i < n; ++i)
  {
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(idx[i] - 1) = 1.0;
    // lower endpoint
    double lower_ub = par0[idx[i] - 1];
    double lower_lb = par0[idx[i] - 1] - 1.0 / std::log(x.rows());
    // lower bound for lower endpoint
    while (2.0 * CEL(method, par0, x, lhs,
                     Eigen::Matrix<double, 1, 1>(lower_lb), maxit, maxit_l, tol,
                     tol_l, gamma, test_th, w)
                     .nllr <=
           cutoff)
    {
      lower_ub = lower_lb;
      lower_lb -= 1.0 / std::log(x.rows());
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > tol)
    {
      if (2.0 * CEL(method, par0, x, lhs,
                    Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2.0),
                    maxit, maxit_l, tol, tol_l, gamma, test_th, w)
                    .nllr >
          cutoff)
      {
        lower_lb = (lower_lb + lower_ub) / 2.0;
      }
      else
      {
        lower_ub = (lower_lb + lower_ub) / 2.0;
      }
    }
    ci(i, 0) = lower_ub;

    // upper endpoint
    double upper_lb = par0[idx[i] - 1];
    double upper_ub = par0[idx[i] - 1] + 1.0 / std::log(x.rows());
    // upper bound for upper endpoint
    while (2.0 * CEL(method, par0, x, lhs,
                     Eigen::Matrix<double, 1, 1>(upper_ub), maxit, maxit_l, tol,
                     tol_l, gamma, test_th, w)
                     .nllr <=
           cutoff)
    {
      upper_lb = upper_ub;
      upper_ub += 1.0 / std::log(x.rows());
    }
    // approximate upper bound by numerical search
    while (upper_ub - upper_lb > tol)
    {
      if (2.0 * CEL(method, par0, x, lhs,
                    Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2.0),
                    maxit, maxit_l, tol, tol_l, gamma, test_th, w)
                    .nllr >
          cutoff)
      {
        upper_ub = (upper_lb + upper_ub) / 2.0;
      }
      else
      {
        upper_lb = (upper_lb + upper_ub) / 2.0;
      }
    }
    ci(i, 1) = upper_lb;
  }
  return ci;
}
