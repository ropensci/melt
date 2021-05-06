#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
Rcpp::NumericVector plog(Rcpp::NumericVector x, double threshold);
Rcpp::NumericVector dplog(Rcpp::NumericVector x, double threshold);
Rcpp::NumericVector d2plog(Rcpp::NumericVector x, double threshold);
double dfp1dcpp(arma::vec gr);
void theta2lambda_bibd(const arma::rowvec &theta,
                       arma::vec &lambda, bool &convex_hull_check,
                       const arma::mat &x, const arma::mat &c,
                       const int &n, const int &p,
                       int maxit = 100, double abstol = 1e-8);
void lambda2theta_bibd(arma::rowvec &theta, arma::vec lambda,
                       const arma::mat &x, const arma::mat &c,
                       const std::vector<int> &pair, double gamma,
                       const int &n, const int &p,
                       int maxit, double abstol);
std::vector<std::vector<int>> all_pairs(const int &p);
arma::mat cov_estimator(const arma::mat &x,
                        const arma::mat &c);
arma::mat g_mean(const arma::vec &theta, arma::mat x);
arma::mat g_ibd(const arma::vec& theta, const arma::mat& x, const arma::mat& c);
arma::vec linear_projection(const arma::vec &theta,
                            const arma::mat &L,
                            const arma::vec &rhs);
arma::vec lambda2theta_ibd(const arma::vec &lambda,
                           const arma::vec &theta,
                           const arma::mat &g,
                           const arma::mat &c,
                           const double &gamma);
#endif
