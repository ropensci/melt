#ifndef EL_NULL_TRANSFORM_H_
#define EL_NULL_TRANSFORM_H_

#include <RcppEigen.h>
#include <functional>

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &)>
  transform_x_fn(const std::string method);

Eigen::MatrixXd x0_mean(const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::VectorXd> &par,
                        const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_sd(const Eigen::Ref<const Eigen::MatrixXd> &x,
                      const Eigen::Ref<const Eigen::VectorXd> &par,
                      const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_lm(const Eigen::Ref<const Eigen::MatrixXd> &x,
                      const Eigen::Ref<const Eigen::VectorXd> &par,
                      const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_gauss_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par,
                             const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_gauss_inverse(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                 const Eigen::Ref<const Eigen::VectorXd> &par,
                                 const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_bin_logit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par,
                             const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_bin_probit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                              const Eigen::Ref<const Eigen::VectorXd> &par,
                              const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_bin_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par,
                           const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_poi_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par,
                           const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_poi_identity(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                const Eigen::Ref<const Eigen::VectorXd> &par,
                                const Eigen::Ref<const Eigen::VectorXd> &est);

Eigen::MatrixXd x0_poi_sqrt(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par,
                            const Eigen::Ref<const Eigen::VectorXd> &est);
#endif
