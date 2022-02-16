#ifndef EL_METHODS_H_
#define EL_METHODS_H_

#include "EL.h"

Eigen::MatrixXd EL_confint(const Eigen::Map<Eigen::MatrixXd>& x,
                           const std::string type,
                           const Eigen::Map<Eigen::VectorXd>& init,
                           const double cutoff,
                           const std::vector<int>& idx,
                           const int maxit,
                           const double abstol,
                           const Rcpp::Nullable<double> threshold);

Eigen::MatrixXd EL_confint2(const Eigen::Map<Eigen::MatrixXd>& x,
                            const std::string type,
                            const Eigen::VectorXd par0,
                            const double cutoff,
                            const int maxit,
                            const double abstol,
                            const Rcpp::Nullable<double> threshold);
#endif
