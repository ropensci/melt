#include "apply_null_transformation.h"
#include "helpers.h"
#include "utils.h"
#include "EL.h"
#include <RcppEigen.h>
#include <boost/random/uniform_int_distribution.hpp>
#include <dqrng.h>
#include <xoshiro.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector compute_bootstrap_calibration(
    const double alpha,
    const double statistic,
    const int B,
    const int seed,
    const int nthreads,
    const std::string method,
    const Eigen::Map<Eigen::MatrixXd> &x,
    const Eigen::Map<Eigen::VectorXd> &par,
    const Eigen::Map<Eigen::VectorXd> &est,
    const int maxit_l,
    const double tol_l,
    const Rcpp::Nullable<double> th,
    const Eigen::Map<Eigen::ArrayXd> &w)
{
  // estimating function
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd> &,
      const Eigen::Ref<const Eigen::VectorXd> &)>
      g_fn = set_g_fn(method);
  // null transformation
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd> &,
      const Eigen::Ref<const Eigen::VectorXd> &,
      const Eigen::Ref<const Eigen::VectorXd> &)>
      x0_fn = transform_x_fn(method);
  const Eigen::MatrixXd g0 = g_fn(x0_fn(x, par, est), par);
  // initialize seed
  dqrng::xoshiro256plus gen(seed);
  // discrete uniform distribution
  boost::random::uniform_int_distribution<> u(0, g0.rows() - 1);
  // bootstrap statistics
  std::vector<double> boot_statistic(B);
  const double test_th = set_threshold(par.size(), th);
  #ifdef _OPENMP
  #pragma omp parallel num_threads(nthreads)
  {
  #endif
    // make thread local copy of gen
    dqrng::xoshiro256plus lgen(gen);
    #ifdef _OPENMP
    // advance gen by 1 ... nthreads jumps
    lgen.jump(omp_get_thread_num() + 1);
    #pragma omp for
    #endif
    for (int i = 0; i < B; ++i)
    {
      Eigen::MatrixXd boot_g(g0.rows(), g0.cols());
      for (int j = 0; j < g0.rows(); ++j)
      {
        boot_g.row(j) = g0.row(u(lgen));
      }
      const EL el(boot_g, maxit_l, tol_l, test_th, w);
      boot_statistic[i] = 2.0 * el.nllr;
    }
  #ifdef _OPENMP
  }
  #endif
  // p-value
  boot_statistic.erase(
      std::remove_if(std::begin(boot_statistic), std::end(boot_statistic),
                     [](const double &x) { return std::isnan(x); }),
      std::end(boot_statistic));
  const double pval =
      count_if(boot_statistic.begin(), boot_statistic.end(),
               [statistic](const double &x) { return (x > statistic); }) /
      static_cast<double>(boot_statistic.size());
  const Rcpp::NumericVector out =
      Rcpp::NumericVector::create(Rcpp::Named("cv") = compute_quantile(
                                      Rcpp::wrap(boot_statistic), 1 - alpha),
                                  Rcpp::Named("pval") = pval);
  return out;
}
