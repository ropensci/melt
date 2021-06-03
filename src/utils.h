#ifndef EL_UTILS_H_
#define EL_UTILS_H_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
  double p_value;
  int iterations;
  bool convergence;
};

arma::vec plog(const arma::vec& x, const double threshold);
arma::vec plog_old(const arma::vec& x, const double threshold);
arma::vec dplog(const arma::vec& x, const double threshold);
arma::vec dplog_old(const arma::vec& x, const double threshold);
arma::vec d2plog(const arma::vec& x, const double threshold);
arma::vec d2plog_old(const arma::vec& x, const double threshold);
EL getEL(const arma::mat& g,
         const int maxit = 100,
         const double abstol = 1e-8);
std::vector<std::array<int, 2>> all_pairs(const int p);
arma::mat g_mean(const arma::vec& theta, arma::mat x);
arma::vec linear_projection(const arma::vec& theta,
                            const arma::mat& L,
                            const arma::vec& rhs);
arma::mat bootstrap_sample(const arma::mat& x);
#endif
