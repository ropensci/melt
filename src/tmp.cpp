#include <xoshiro.h>
#include <dqrng.h>
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <random>
#include <vector>

// [[Rcpp::export]]
double tt2(const int B, const int seed, const int nthreads,
           const Eigen::MatrixXd mat) {
  const int n = mat.rows();
  const int p = mat.cols();
  // initialize seed
  dqrng::xoshiro256plus gen(seed);

  std::uniform_int_distribution<> dist(0, n - 1);

  std::vector<double> bootstrap2(B);
  #ifdef _OPENMP
  #pragma omp parallel num_threads(nthreads)
  {
  #endif
    dqrng::xoshiro256plus lgen(gen);      // make thread local copy of rng
    #ifdef _OPENMP
    lgen.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... ncores jumps
    #pragma omp for
    #endif
    for (int i = 0; i < B; ++i) {
      // place holder for the bootstrap data matrix
      Eigen::MatrixXd out2(n, p);
      for (int j = 0; j < n; ++j) {
        out2.row(j) = mat.row(dist(lgen));
        }
      // once out2 matrix is filled, compute EL with it and fill bootstrap2
      bootstrap2[i] = out2.mean();
    }
  #ifdef _OPENMP
  }
  #endif

  Rcpp::Environment stats("package:stats");
  Rcpp::Function quantile = stats["quantile"];
  return Rcpp::as<double>(quantile(bootstrap2, Rcpp::Named("probs") = 0.95));
}