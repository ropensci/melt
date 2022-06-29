#ifndef EL_HELPERS_H_
#define EL_HELPERS_H_

#include <RcppEigen.h>

double setThreshold(const int p, const Rcpp::Nullable<double> th);
double setStep(const int n, const Rcpp::Nullable<double> step);

#endif
