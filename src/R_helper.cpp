#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
int get_max_threads() {
  #ifdef _OPENMP
    return omp_get_max_threads();
  #else
    return 1;
  #endif
};

// [[Rcpp::export]]
int get_rank(const Eigen::Map<Eigen::MatrixXd>& x) {
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  return lu_decomp.rank();
}
