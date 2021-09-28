#ifndef EL_UTILS_gbd_H_
#define EL_UTILS_gbd_H_

#include "EL.h"
#include "utils.h"
#include <omp.h>

#include <progress.hpp>
#include "nb_progress_bar.h"



Eigen::MatrixXd g_gbd(const Eigen::Ref<const Eigen::VectorXd>& theta,
                      const Eigen::Ref<const Eigen::MatrixXd>& x,
                      const Eigen::Ref<const Eigen::MatrixXd>& c);

Eigen::MatrixXd cov_gbd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                        const Eigen::Ref<const Eigen::MatrixXd>& c);

Eigen::VectorXd lambda2theta_gbd(const Eigen::Ref<const Eigen::VectorXd>& lambda,
                                 const Eigen::Ref<const Eigen::VectorXd>& theta,
                                 const Eigen::Ref<const Eigen::MatrixXd>& g,
                                 const Eigen::Ref<const Eigen::MatrixXd>& c,
                                 const double gamma);

void lambda2theta_void(
        const Eigen::Ref<const Eigen::VectorXd>& lambda,
        Eigen::Ref<Eigen::VectorXd> theta,
        const Eigen::Ref<const Eigen::MatrixXd>& g,
        const Eigen::Ref<const Eigen::MatrixXd>& c,
        const double gamma);

Eigen::VectorXd approx_lambda_gbd(
        const Eigen::Ref<const Eigen::MatrixXd>& g0,
        const Eigen::Ref<const Eigen::MatrixXd>& c,
        const Eigen::Ref<const Eigen::VectorXd>& theta0,
        const Eigen::Ref<const Eigen::VectorXd>& theta1,
        const Eigen::Ref<const Eigen::VectorXd>& lambda0);

std::array<double, 2> pair_confidence_interval_gbd(
    const Eigen::Ref<const Eigen::VectorXd>& theta0,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const double threshold,
    const double init,
    const double cutoff);

Eigen::MatrixXd rmvn(const Eigen::MatrixXd& x, const int n);

Eigen::ArrayXd bootstrap_statistics_pairwise_AMC(
        const Eigen::Ref<const Eigen::MatrixXd>& x,
        const Eigen::Ref<const Eigen::MatrixXd>& c,
        const int k,
        const std::vector<std::array<int, 2>>& pairs,
        const int B,
        const double level);

Eigen::ArrayXd bootstrap_statistics_pairwise_NB(
        const Eigen::Ref<const Eigen::MatrixXd>& x,
        const Eigen::Ref<const Eigen::MatrixXd>& c,
        const int k,
        const std::vector<std::array<int, 2>>& pairs,
        const int B,
        const double level,
        const int nthread,
        const bool progress,
        const double threshold,
        const int maxit,
        const double abstol);

minEL test_gbd_EL(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                  const Eigen::Ref<const Eigen::MatrixXd>& x,
                  const Eigen::Ref<const Eigen::MatrixXd>& c,
                  const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                  const Eigen::Ref<const Eigen::VectorXd>& rhs,
                  const double threshold,
                  const int maxit = 1000,
                  const double abstol = 1e-8);

double test_nlogLR(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                   const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const double threshold,
                   const int maxit = 1000,
                   const double abstol = 1e-8);

double test_nlogLR(const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const double threshold,
                   const int maxit = 1000,
                   const double abstol = 1e-8);
#endif
