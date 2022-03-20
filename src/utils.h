#ifndef EL_UTILS_H_
#define EL_UTILS_H_

#include <RcppEigen.h>

double th_nlogLR(const int p, const Rcpp::Nullable<double> threshold);

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::VectorXd>& par,
                       const Eigen::Ref<const Eigen::MatrixXd>& x);

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::VectorXd>& par,
                     const Eigen::Ref<const Eigen::MatrixXd>& data);

Eigen::VectorXd gr_nlogLR_lm(
        const Eigen::Ref<const Eigen::VectorXd>& lambda,
        const Eigen::Ref<const Eigen::MatrixXd>& g,
        const Eigen::Ref<const Eigen::MatrixXd>& data);
#endif
