#include "methods.h"

// [[Rcpp::export]]
Eigen::MatrixXd EL_confint(const Eigen::Map<Eigen::MatrixXd>& x,
                           const std::string type,
                           const Eigen::Map<Eigen::VectorXd>& init,
                           const double cutoff,
                           const int maxit,
                           const double abstol) {
  const int p = x.cols();
  Eigen::MatrixXd ci(p, 2);
  for (int j = 0; j < p; ++j) {
    // upper endpoint
    double upper_lb = init[j];
    double upper_size = 1;
    double upper_ub = init[j] + upper_size;
    // upper bound for upper endpoint
    while (2 * EL2(Eigen::Matrix<double, 1, 1>(upper_ub), x, type, maxit,
                   abstol).nlogLR <= cutoff) {
      upper_lb = upper_ub;
      upper_ub += upper_size;
    }
    // approximate upper bound by numerical search
    while (upper_ub - upper_lb > 1e-04) {
      if (2 * EL2(Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2), x,
                  type, maxit, abstol).nlogLR > cutoff) {
        upper_ub = (upper_lb + upper_ub) / 2;
      } else {
        upper_lb = (upper_lb + upper_ub) / 2;
      }
    }
    // lower endpoint
    double lower_ub = init[j];
    double lower_size = 1;
    double lower_lb = init[j] - lower_size;

    // lower bound for lower endpoint
    while (2 * EL2(Eigen::Matrix<double, 1, 1>(lower_lb), x, type, maxit,
                   abstol).nlogLR <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= lower_size;
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > 1e-04) {
      if (2 * EL2(Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2), x,
                  type, maxit, abstol).nlogLR > cutoff) {
        lower_lb = (lower_lb + lower_ub) / 2;
      } else {
        lower_ub = (lower_lb + lower_ub) / 2;
      }
    }
    ci(j, 0) = lower_lb;
    ci(j, 1) = upper_ub;
  }
  return ci;
}

