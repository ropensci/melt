#ifndef EL_UTILS_IBD_H_
#define EL_UTILS_IBD_H_

#include "utils.h"
#include <RcppArmadillo.h>
#include <omp.h>

arma::mat g_ibd(const arma::vec& theta, const arma::mat& x, const arma::mat& c);
arma::mat cov_ibd(const arma::mat& x,
                  const arma::mat& c,
                  const bool adjust);
arma::vec lambda2theta_ibd(const arma::vec& lambda,
                           const arma::vec& theta,
                           const arma::mat& g,
                           const arma::mat& c,
                           const double gamma);
arma::vec approx_lambda_ibd(const arma::mat& x,
                            const arma::mat& c,
                            const arma::vec& theta0,
                            const arma::vec& theta1,
                            const arma::vec& lambda0);
arma::mat centering_ibd(arma::mat x);

arma::umat block_bootstrap_index(const arma::mat& x,
                                 const arma::mat& c,
                                 const int B,
                                 const bool approx_lambda,
                                 const int maxit,
                                 const double abstol);

double cutoff_pairwise_PB_ibd(const arma::mat& x,
                              const arma::mat& c,
                              const int B,
                              const double level,
                              const bool adjust);
double cutoff_pairwise_NPB_ibd(const arma::mat& x,
                               const arma::mat& c,
                               const int B,
                               const double level,
                               const bool block_bootstrap,
                               const bool approx_lambda,
                               const int ncores,
                               const int maxit,
                               const double abstol);
std::array<double, 2> pair_confidence_interval_ibd(const arma::mat& x,
                                                   const arma::mat& c,
                                                   const arma::mat& L,
                                                   const bool approx_lambda,
                                                   const double init,
                                                   const double threshold);
minEL test_ibd_EL(const arma::mat& x,
                  const arma::mat& c,
                  const arma::mat& L,
                  const arma::vec& rhs,
                  const bool approx_lambda,
                  const int maxit = 1000,
                  const double abstol = 1e-8);
#endif
