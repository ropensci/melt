#ifndef EL_UTILS_H_
#define EL_UTILS_H_

#include <RcppEigen.h>

int get_max_threads();

int get_rank(const Eigen::Map<Eigen::MatrixXd> &x);

double set_threshold(const int p, const Rcpp::Nullable<double> th);

double set_step(const int n, const Rcpp::Nullable<double> step);

double compute_quantile(const Rcpp::NumericVector &x, const double prob);
#endif
