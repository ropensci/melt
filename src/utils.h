#ifndef EL_UTILS_H_
#define EL_UTILS_H_

#include <RcppEigen.h>

double step_nloglr(const int n, const Rcpp::Nullable<double> step);
double th_nloglr(const int p, const Rcpp::Nullable<double> th);

Eigen::ArrayXd inverse_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x);
Eigen::ArrayXd log_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x);
Eigen::ArrayXd logit_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x);
Eigen::ArrayXd probit_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x);


Eigen::VectorXd mele_mean(const Eigen::Ref<const Eigen::MatrixXd>& x,
                          const Eigen::Ref<const Eigen::ArrayXd>& w);
Eigen::VectorXd mele_lm(const Eigen::Ref<const Eigen::MatrixXd>& data,
                        const Eigen::Ref<const Eigen::ArrayXd>& w);


Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::MatrixXd>& x,
                       const Eigen::Ref<const Eigen::VectorXd>& par);
Eigen::VectorXd gr_nloglr_mean(const Eigen::Ref<const Eigen::VectorXd>& l,
                               const Eigen::Ref<const Eigen::MatrixXd>& g,
                               const Eigen::Ref<const Eigen::MatrixXd>& x,
                               const Eigen::Ref<const Eigen::VectorXd>& par,
                               const Eigen::Ref<const Eigen::ArrayXd>& w,
                               const bool weighted);


// Gaussian family
Eigen::MatrixXd g_gauss_log(const Eigen::Ref<const Eigen::MatrixXd>& data,
                               const Eigen::Ref<const Eigen::VectorXd>& par);
Eigen::VectorXd gr_nloglr_gauss_log(
                const Eigen::Ref<const Eigen::VectorXd>& l,
                const Eigen::Ref<const Eigen::MatrixXd>& g,
                const Eigen::Ref<const Eigen::MatrixXd>& data,
                const Eigen::Ref<const Eigen::VectorXd>& par,
                const Eigen::Ref<const Eigen::ArrayXd>& w,
                const bool weighted);

Eigen::MatrixXd g_gauss_inverse(const Eigen::Ref<const Eigen::MatrixXd>& data,
                                const Eigen::Ref<const Eigen::VectorXd>& par);
Eigen::VectorXd gr_nloglr_gauss_inverse(
                const Eigen::Ref<const Eigen::VectorXd>& l,
                const Eigen::Ref<const Eigen::MatrixXd>& g,
                const Eigen::Ref<const Eigen::MatrixXd>& data,
                const Eigen::Ref<const Eigen::VectorXd>& par,
                const Eigen::Ref<const Eigen::ArrayXd>& w,
                const bool weighted);


// quantile function
double quantileRcpp(const Rcpp::NumericVector& x, double prob);

#endif
