#ifndef EL_UTILS_H_
#define EL_UTILS_H_

#include <RcppEigen.h>

double th_nloglr(const int p, const Rcpp::Nullable<double> th);

Eigen::VectorXd mele_mean(const Eigen::Ref<const Eigen::MatrixXd>& x,
                          const Eigen::Ref<const Eigen::ArrayXd>& w);
Eigen::VectorXd mele_lm(const Eigen::Ref<const Eigen::MatrixXd>& data,
                        const Eigen::Ref<const Eigen::ArrayXd>& w);

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::MatrixXd>& x,
                       const Eigen::Ref<const Eigen::VectorXd>& par);

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::MatrixXd>& x,
                     const Eigen::Ref<const Eigen::VectorXd>& par);

Eigen::VectorXd gr_nloglr_mean(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w);

Eigen::VectorXd wgr_nloglr_mean(
        const Eigen::Ref<const Eigen::VectorXd>& l,
        const Eigen::Ref<const Eigen::MatrixXd>& g,
        const Eigen::Ref<const Eigen::MatrixXd>& data,
        const Eigen::Ref<const Eigen::ArrayXd>& w);

Eigen::VectorXd gr_nloglr_lm(
        const Eigen::Ref<const Eigen::VectorXd>& l,
        const Eigen::Ref<const Eigen::MatrixXd>& g,
        const Eigen::Ref<const Eigen::MatrixXd>& data,
        const Eigen::Ref<const Eigen::VectorXd>& par,
        const Eigen::Ref<const Eigen::ArrayXd>& w);

Eigen::VectorXd wgr_nloglr_lm(
        const Eigen::Ref<const Eigen::VectorXd>& l,
        const Eigen::Ref<const Eigen::MatrixXd>& g,
        const Eigen::Ref<const Eigen::MatrixXd>& data,
        const Eigen::Ref<const Eigen::ArrayXd>& w);








Eigen::ArrayXd logit_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x);

Eigen::MatrixXd g_logit(const Eigen::Ref<const Eigen::MatrixXd>& data,
                        const Eigen::Ref<const Eigen::VectorXd>& par);

Eigen::VectorXd gr_nloglr_logit(
        const Eigen::Ref<const Eigen::VectorXd>& l,
        const Eigen::Ref<const Eigen::MatrixXd>& g,
        const Eigen::Ref<const Eigen::MatrixXd>& data,
        const Eigen::Ref<const Eigen::VectorXd>& par,
        const Eigen::Ref<const Eigen::ArrayXd>& w);



Eigen::ArrayXd log_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x);
Eigen::MatrixXd g_gaussian_log(const Eigen::Ref<const Eigen::MatrixXd>& data,
                               const Eigen::Ref<const Eigen::VectorXd>& par);
Eigen::VectorXd gr_nloglr_gaussian_log(
        const Eigen::Ref<const Eigen::VectorXd>& l,
        const Eigen::Ref<const Eigen::MatrixXd>& g,
        const Eigen::Ref<const Eigen::MatrixXd>& data,
        const Eigen::Ref<const Eigen::VectorXd>& par);










#endif
