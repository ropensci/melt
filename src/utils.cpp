#include <RcppArmadillo.h>
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::NumericVector plogcpp(Rcpp::NumericVector x, double threshold) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = std::log(x[i]);
    } else {
      out[i] = std::log(threshold) - 1.5 + 2 * std::pow(threshold, -1) * x[i] -
        std::pow(x[i] / threshold, 2) / 2;
    }
  }
  return out;
}

Rcpp::NumericVector dplogcpp(Rcpp::NumericVector x, double threshold) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = std::pow(x[i], -1);
    } else {
      out[i] = 2 * std::pow(threshold, -1) - x[i] * std::pow(threshold, -2);
    }
  }
  return out;
}

Rcpp::NumericVector d2plogcpp(Rcpp::NumericVector x, double threshold) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = -std::pow(x[i], -2);
    } else {
      out[i] = -std::pow(threshold, -2);
    }
  }
  return out;
}

double dfp1dcpp(Rcpp::NumericVector gr) {
  /* This function can be modified to for a general direction finding problem
     in higher dimensions
  */
  /* Prototype:
  arma::mat LHS = {{1, 0, 1}, {0, 1, -1}, {1, -1, 0}};
  arma::vec RHS = {0, 0, 0};
  RHS.subvec(0, 1) = -gr;
  arma::vec sol = solve(LHS, RHS);
  return sol(0);
  */
  return -sum(gr) / 2;
}
