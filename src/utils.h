#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
Rcpp::NumericVector plogcpp(Rcpp::NumericVector x, double threshold);
Rcpp::NumericVector dplogcpp(Rcpp::NumericVector x, double threshold);
Rcpp::NumericVector d2plogcpp(Rcpp::NumericVector x, double threshold);
double dfp1dcpp(arma::vec gr);

#endif
