#ifndef EL_UTILS_IBD_H
#define EL_UTILS_IBD_H

#include "utils.h"
#include <RcppArmadillo.h>
arma::mat g_ibd(const arma::vec& theta, const arma::mat& x, const arma::mat& c);
arma::mat cov_ibd(const arma::mat& x, const arma::mat& c, const bool adjust);
arma::vec lambda2theta_ibd(const arma::vec& lambda,
                           const arma::vec& theta,
                           const arma::mat& g,
                           const arma::mat& c,
                           const double gamma);
double cutoff_pairwise_PB_ibd(const arma::mat& x,
                              const arma::mat& c,
                              const int B,
                              const double level,
                              const bool approx_lambda,
                              const bool adjust);
arma::mat centering_ibd(arma::mat x);
minEL test_ibd_EL(const arma::mat& x,
                  const arma::mat& c,
                  const arma::mat& L,
                  const arma::vec& rhs,
                  const bool approx_lambda,
                  const int maxit = 1000,
                  const double abstol = 1e-8);
std::vector<double> pair_confidence_interval_ibd(const arma::mat& x,
                                                 const arma::mat& c,
                                                 const arma::mat& L,
                                                 const bool approx_lambda,
                                                 const double init,
                                                 const double threshold);
double cutoff_pairwise_NPB_ibd(const arma::mat& x,
                               const int B,
                               const double level,
                               const bool approx_lambda,
                               const int maxit,
                               const double abstol);
double quantile_pairwise_NPB_ibd(const arma::mat& x,
                                 const int B,
                                 const double level,
                                 const int maxit,
                                 const double abstol);
#endif
