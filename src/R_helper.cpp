#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
int max_threads_()
{
  int n = 1;
  #ifdef _OPENMP
  n = omp_get_max_threads();
  #endif
  return n;
};

// [[Rcpp::export]]
int get_rank_(const Eigen::Map<Eigen::MatrixXd>& x)
{
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  return lu_decomp.rank();
}
