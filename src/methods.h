#ifndef EL_METHODS_H_
#define EL_METHODS_H_

#include "EL.h"
Eigen::MatrixXd EL_confint(const Eigen::Map<Eigen::MatrixXd>& x,
                               const std::string type,
                               const Eigen::Map<Eigen::VectorXd>& init,
                               const double cutoff,
                               const std::vector<int>& idx,
                               const int maxit,
                               const double abstol);
#endif
