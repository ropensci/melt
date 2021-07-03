#ifndef EL_UTILS_IBD_H_
#define EL_UTILS_IBD_H_

#include "utils.h"
#include <omp.h>

Eigen::MatrixXd g_ibd(const Eigen::Ref<const Eigen::VectorXd>& theta,
                      const Eigen::Ref<const Eigen::MatrixXd>& x,
                      const Eigen::Ref<const Eigen::MatrixXd>& c);

Eigen::MatrixXd cov_ibd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                        const Eigen::Ref<const Eigen::MatrixXd>& c);

Eigen::VectorXd lambda2theta_ibd(const Eigen::Ref<const Eigen::VectorXd>& lambda,
                                 const Eigen::Ref<const Eigen::VectorXd>& theta,
                                 const Eigen::Ref<const Eigen::MatrixXd>& g,
                                 const Eigen::Ref<const Eigen::MatrixXd>& c,
                                 const double gamma);

Eigen::VectorXd approx_lambda_ibd(
        const Eigen::Ref<const Eigen::MatrixXd>& g0,
        const Eigen::Ref<const Eigen::MatrixXd>& c,
        const Eigen::Ref<const Eigen::VectorXd>& theta0,
        const Eigen::Ref<const Eigen::VectorXd>& theta1,
        const Eigen::Ref<const Eigen::VectorXd>& lambda0);

std::array<double, 2> pair_confidence_interval_ibd(
    const Eigen::Ref<const Eigen::VectorXd>& theta0,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const bool approx_lambda,
    const double init,
    const double threshold);

Eigen::MatrixXd centering_ibd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                              const Eigen::Ref<const Eigen::MatrixXd>& c);

Eigen::MatrixXd rmvn(const Eigen::MatrixXd& x, const int n);

double cutoff_pairwise_PB_ibd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                              const Eigen::Ref<const Eigen::MatrixXd>& c,
                              const std::vector<std::array<int, 2>>& pairs,
                              const int B,
                              const double level);
double cutoff_pairwise_NB_ibd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                              const Eigen::Ref<const Eigen::MatrixXd>& c,
                              const int B,
                              const double level,
                              const bool block_bootstrap,
                              const bool approx_lambda,
                              const int ncores,
                              const int maxit,
                              const double abstol);

// initial value given
minEL test_ibd_EL(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                  const Eigen::Ref<const Eigen::MatrixXd>& x,
                  const Eigen::Ref<const Eigen::MatrixXd>& c,
                  const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                  const Eigen::Ref<const Eigen::VectorXd>& rhs,
                  const bool approx_lambda,
                  const int maxit = 1000,
                  const double abstol = 1e-8);
// initial value not given
minEL test_ibd_EL(const Eigen::Ref<const Eigen::MatrixXd>& x,
                  const Eigen::Ref<const Eigen::MatrixXd>& c,
                  const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                  const Eigen::Ref<const Eigen::VectorXd>& rhs,
                  const bool approx_lambda,
                  const int maxit = 1000,
                  const double abstol = 1e-8);
#endif
