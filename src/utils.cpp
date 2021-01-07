#include "utils.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

NumericVector plogcpp(NumericVector x, double threshold) {
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = log(x[i]);
    } else {
      out[i] = log(threshold) - 1.5 + 2 * pow(threshold, -1) * x[i] -
        pow(x[i] / threshold, 2) / 2;
    }
  }
  return out;
}

NumericVector dplogcpp(NumericVector x, double threshold) {
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = pow(x[i], -1);
    } else {
      out[i] = 2 * pow(threshold, -1) - x[i] * pow(threshold, -2);
    }
  }
  return out;
}

NumericVector d2plogcpp(NumericVector x, double threshold) {
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = -pow(x[i], -2);
    } else {
      out[i] = -pow(threshold, -2);
    }
  }
  return out;
}

double dfp1dcpp(vec gr) {
  mat LHS = {{1, 0, 1}, {0, 1, -1}, {1, -1, 0}};
  vec RHS = zeros<vec>(3);
  RHS.subvec(0, 1) = -gr;
  vec sol = solve(LHS, RHS);
  return sol(0);
}
