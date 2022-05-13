#ifndef MHT_UTILS_H_
#define MHT_UTILS_H_

#include "utils.h"

Eigen::RowVectorXd rmvn(const Eigen::Ref<const Eigen::MatrixXd>& sqrt);

Eigen::MatrixXd w_mean(const Eigen::Ref<const Eigen::MatrixXd>& x);

Eigen::MatrixXd dg0_inv(const std::string method,
                        const Eigen::Ref<const Eigen::MatrixXd>& x);

Eigen::MatrixXd cov(const std::string method,
                    const Eigen::Ref<const Eigen::VectorXd>& est,
                    const Eigen::Ref<const Eigen::MatrixXd>& x);

Eigen::MatrixXd ahat(const Eigen::Ref<const Eigen::MatrixXd>& j,
                     const Eigen::Ref<const Eigen::MatrixXd>& w,
                     const Eigen::Ref<const Eigen::MatrixXd>& s);

Eigen::MatrixXd cov2(const std::string method,
                     const Eigen::Map<Eigen::VectorXd>& est,
                     const Eigen::Map<Eigen::MatrixXd>& x);

#endif
