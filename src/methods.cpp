#include "methods.h"

// [[Rcpp::export]]
Eigen::MatrixXd EL_confint(const Eigen::Map<Eigen::MatrixXd>& x,
                           const std::string type,
                           const Eigen::Map<Eigen::VectorXd>& init,
                           const double cutoff,
                           const std::vector<int>& idx,
                           const int maxit,
                           const double abstol,
                           const Rcpp::Nullable<double> threshold) {
  const int p = idx.size();
  std::vector<double> ci_vec;
  ci_vec.reserve(2 * p);
  // for (int j = 0; j < p; ++j) {
  for (int j : idx) {
    // lower endpoint
    double lower_ub = init[j - 1];
    double lower_size = 1.0;
    double lower_lb = init[j - 1] - lower_size;
    // lower bound for lower endpoint
    while (2.0 * EL2(Eigen::Matrix<double, 1, 1>(lower_lb), x, type,
                     maxit, abstol, th_nlogLR(1, threshold)).nlogLR <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= lower_size;
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > 1e-04) {
      if (2.0 * EL2(Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2), x,
                    type, maxit, abstol, th_nlogLR(1, threshold)).nlogLR > cutoff) {
        lower_lb = (lower_lb + lower_ub) / 2;
      } else {
        lower_ub = (lower_lb + lower_ub) / 2;
      }
    }
    ci_vec.push_back(lower_lb);

    // upper endpoint
    double upper_lb = init[j - 1];
    double upper_size = 1.0;
    double upper_ub = init[j - 1] + upper_size;
    // upper bound for upper endpoint
    while (2.0 * EL2(Eigen::Matrix<double, 1, 1>(upper_ub), x, type,
                     maxit, abstol, th_nlogLR(1, threshold)).nlogLR <= cutoff) {
      upper_lb = upper_ub;
      upper_ub += upper_size;
    }
    // approximate upper bound by numerical search
    while (upper_ub - upper_lb > 1e-04) {
      if (2.0 * EL2(Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2), x,
                  type, maxit, abstol, th_nlogLR(1, threshold)).nlogLR > cutoff) {
        upper_ub = (upper_lb + upper_ub) / 2;
      } else {
        upper_lb = (upper_lb + upper_ub) / 2;
      }
    }
    ci_vec.push_back(upper_ub);
  }
  Eigen::MatrixXd ci =
    Eigen::Map<Eigen::MatrixXd>(ci_vec.data(), 2, p);
  return ci.transpose();
}

// [[Rcpp::export]]
Eigen::MatrixXd EL_confint2(const Eigen::Map<Eigen::MatrixXd>& x,
                            const std::string type,
                            const Eigen::VectorXd par0,
                            const double cutoff,
                            const int maxit,
                            const double abstol,
                            const Rcpp::Nullable<double> threshold) {
  const int p = par0.size();
  std::vector<double> ci_vec;
  ci_vec.reserve(2 * p);
  for (int j = 0; j < 2; ++j) {
  // for (int j : idx) {
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(j) = 1.0;

    // lower endpoint
    double lower_ub = par0[j];
    double lower_size = 1.0;
    double lower_lb = par0[j] - lower_size;
    Rcpp::Rcout << lower_lb << "\n";
    Rcpp::Rcout << lower_ub << "\n";
    // lower bound for lower endpoint
    while (2.0 * EL2(par0, x, type, lhs, Eigen::Matrix<double, 1, 1>(lower_lb),
                     maxit, abstol, 20000).nlogLR <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= lower_size;
    }
    Rcpp::Rcout << lower_lb << "\n";
    Rcpp::Rcout << lower_ub << "\n";
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > 1e-04) {
      if (2.0 * EL2(par0, x, type, lhs,
                    Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2),
                    maxit, abstol, 20000).nlogLR > cutoff) {
        lower_lb = (lower_lb + lower_ub) / 2;
      } else {
        lower_ub = (lower_lb + lower_ub) / 2;
      }
    }
    ci_vec.push_back(lower_lb);

    // upper endpoint
    double upper_lb = par0[j];
    double upper_size = 1.0;
    double upper_ub = par0[j] + upper_size;

    // upper bound for upper endpoint
    while (2.0 * EL2(par0, x, type, lhs, Eigen::Matrix<double, 1, 1>(upper_ub),
                     maxit, abstol, th_nlogLR(1, threshold)).nlogLR <= cutoff) {
      upper_lb = upper_ub;
      upper_ub += upper_size;
    }
    // approximate upper bound by numerical search
    while (upper_ub - upper_lb > 1e-04) {
      if (2.0 * EL2(par0, x, type, lhs,
                    Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2),
                    maxit, abstol, th_nlogLR(1, threshold)).nlogLR > cutoff) {
        upper_ub = (upper_lb + upper_ub) / 2;
      } else {
        upper_lb = (upper_lb + upper_ub) / 2;
      }
    }
    ci_vec.push_back(upper_ub);
  }
  Eigen::MatrixXd ci =
    Eigen::Map<Eigen::MatrixXd>(ci_vec.data(), 2, p);
  return ci.transpose();
}

