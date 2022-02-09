#ifndef EL_UTILS_H_
#define EL_UTILS_H_

#include <RcppEigen.h>

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::VectorXd>& par,
                       const Eigen::Ref<const Eigen::MatrixXd>& x);

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::VectorXd>& beta,
                     const Eigen::Ref<const Eigen::MatrixXd>& x,
                     const Eigen::Ref<const Eigen::VectorXd>& y);

Eigen::MatrixXd estimating_func(
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::VectorXd>&,
                                  const Eigen::Ref<const Eigen::MatrixXd>&)> func);

std::vector<std::array<int, 2>> comparison_pairs(
        const int p, const int control);

Eigen::VectorXd linear_projection(
    const Eigen::Ref<const Eigen::VectorXd>& theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs);
void linear_projection_void(
    Eigen::Ref<Eigen::VectorXd> theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs);

Eigen::MatrixXd bootstrap_sample(const Eigen::Ref<const Eigen::MatrixXd>& x,
                                 const Eigen::Ref<const Eigen::ArrayXi>& index);
#endif
