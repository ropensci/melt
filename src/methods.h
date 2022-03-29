#ifndef EL_METHODS_H_
#define EL_METHODS_H_

#include "EL.h"

Eigen::MatrixXd confint_(const std::string method,
                         const Eigen::Map<Eigen::VectorXd>& par0,
                         const Eigen::Map<Eigen::MatrixXd>& x,
                         const double cutoff,
                         const Rcpp::IntegerVector& idx,
                         const int maxit,
                         const double tol,
                         const Rcpp::Nullable<double> threshold);

Eigen::MatrixXd confint_w_(const std::string method,
                           const Eigen::Map<Eigen::VectorXd>& par0,
                           const Eigen::Map<Eigen::MatrixXd>& x,
                           const Eigen::Map<Eigen::ArrayXd>& w,
                           const double cutoff,
                           const Rcpp::IntegerVector& idx,
                           const int maxit,
                           const double tol,
                           const Rcpp::Nullable<double> th);
#endif
