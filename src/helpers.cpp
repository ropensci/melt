#include "helpers.h"
#include <RcppEigen.h>

double setThreshold(const int p, const Rcpp::Nullable<double> th) {
  return (th.isNull())? 200.0 * p : Rcpp::as<double>(th);
}

double setStep(const int n, const Rcpp::Nullable<double> step) {
  return (step.isNull())? static_cast<double>(1.0 / n) : Rcpp::as<double>(step);
}
