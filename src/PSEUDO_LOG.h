#ifndef PSEUDO_LOG_H_
#define PSEUDO_LOG_H_

#include <RcppEigen.h>

class PSEUDO_LOG {
public:
  Eigen::ArrayXd dplog;
  Eigen::ArrayXd sqrt_neg_d2plog;
  double plog_sum;

  // PSEUDO_LOG(const Eigen::Ref<const Eigen::VectorXd>& x);
  PSEUDO_LOG(Eigen::VectorXd&& x);
  static double sum(Eigen::VectorXd&& x);
  static Eigen::ArrayXd dp(Eigen::VectorXd&& x);
};
#endif
