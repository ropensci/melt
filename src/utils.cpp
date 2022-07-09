#include "utils.h"
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
}

// [[Rcpp::export]]
int get_rank(const Eigen::Map<Eigen::MatrixXd> &x) {
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  return lu_decomp.rank();
}

double set_threshold(const int p, const Rcpp::Nullable<double> th) {
  return (th.isNull()) ? 200.0 * p : Rcpp::as<double>(th);
}

double set_step(const int n, const Rcpp::Nullable<double> step) {
  return
  (step.isNull()) ? static_cast<double>(1.0 / n) : Rcpp::as<double>(step);
}

double compute_quantile(const Rcpp::NumericVector &x, const double prob) {
  Rcpp::Environment stats("package:stats");
  Rcpp::Function quantile = stats["quantile"];
  return Rcpp::as<double>(quantile(x, Rcpp::Named("probs") = prob));
}
