#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
int max_threads()
{
  int n = 1;
  #ifdef _OPENMP
  n = omp_get_max_threads();
  #endif
  return n;
};
