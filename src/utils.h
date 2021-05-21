#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
struct EL {
  double nlogLR;
  arma::vec lambda;
  arma::vec gradient;
  int iterations;
  bool convergence;
};
struct minEL {
  arma::vec theta;
  arma::vec lambda;
  double nlogLR;
  int iterations;
  bool convergence;
};

arma::vec plog(const arma::vec& x, double threshold);
Rcpp::NumericVector plog_old(Rcpp::NumericVector x, double threshold);
arma::vec dplog(const arma::vec& x, double threshold);
Rcpp::NumericVector dplog_old(Rcpp::NumericVector x, double threshold);
arma::vec d2plog(const arma::vec& x, double threshold);
Rcpp::NumericVector d2plog_old(Rcpp::NumericVector x, double threshold);
EL getEL(const arma::mat& g,
         const int& maxit = 100,
         const double& abstol = 1e-8);
std::vector<std::vector<int>> all_pairs(const int& p);
arma::mat g_mean(const arma::vec& theta, arma::mat x);
arma::mat g_ibd(const arma::vec& theta, const arma::mat& x, const arma::mat& c);
arma::mat cov_ibd(const arma::mat& x, const arma::mat& c, const bool adjust);
arma::vec linear_projection(const arma::vec& theta,
                            const arma::mat& L,
                            const arma::vec& rhs);
arma::vec lambda2theta_ibd(const arma::vec& lambda,
                           const arma::vec& theta,
                           const arma::mat& g,
                           const arma::mat& c,
                           const double& gamma);
double threshold_pairwise_ibd(const arma::mat& x,
                              const arma::mat& c,
                              const int& B,
                              const double& level,
                              const bool adjust);
#endif
