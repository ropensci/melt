#ifndef EL_UTILS_gbd_H_
#define EL_UTILS_gbd_H_

#include "EL.h"
#include "utils.h"

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

Eigen::MatrixXd rmvn(const Eigen::MatrixXd& x, const int n);

minEL test_gbd_EL(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                  const Eigen::Ref<const Eigen::MatrixXd>& x,
                  const Eigen::Ref<const Eigen::MatrixXd>& c,
                  const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                  const Eigen::Ref<const Eigen::VectorXd>& rhs,
                  const int maxit,
                  const double abstol,
                  const double threshold);

double test_nlogLR(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                   const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const int maxit,
                   const double abstol,
                   const double threshold);

double test_nlogLR(const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const int maxit,
                   const double abstol,
                   const double threshold);
#endif
