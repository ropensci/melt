#include "utils.h"
#include "EL.h"
#include <xoshiro.h>
#include <dqrng.h>
#include <boost/random/uniform_int_distribution.hpp>
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <random>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector boot_(const int B,
                          const int seed,
                          const int nthreads,
                          const std::string method,
                          const Eigen::Map<Eigen::MatrixXd>& x,
                          const Eigen::Map<Eigen::VectorXd>& par,
                          const int maxit_l,
                          const double tol_l,
                          const Rcpp::Nullable<double> th,
                          const Eigen::Map<Eigen::ArrayXd>& wt) {
  // estimating function
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&)> g_fn = g_fn2(method);
  // null transformation
  const Eigen::MatrixXd g0 =
    g_fn(x, par).rowwise() - g_fn(x, par).colwise().mean();
  const int n = g0.rows();
  const int p = g0.cols();
  // initialize seed
  dqrng::xoshiro256plus gen(seed);
  // discrete uniform distribution
  boost::random::uniform_int_distribution<> u(0, n - 1);

  std::vector<double> boot_statistic(B);
  const double test_th = th_nloglr(par.size(), th);
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
    for (int i = 0; i < B; ++i) {
      Eigen::MatrixXd boot_g(n, p);
      for (int j = 0; j < n; ++j) {
        boot_g.row(j) = g0.row(u(lgen));
      }
      const EL el(boot_g, maxit_l, tol_l, test_th, wt);
      boot_statistic[i] = 2.0 * el.nllr;
    }
  #ifdef _OPENMP
  }
  #endif
  const Rcpp::NumericVector out = Rcpp::wrap(boot_statistic);
  return out;
}
