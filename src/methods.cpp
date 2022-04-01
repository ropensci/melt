#include "methods.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
Eigen::MatrixXd confint_(const std::string method,
                         const Eigen::Map<Eigen::VectorXd>& par0,
                         const Eigen::Map<Eigen::MatrixXd>& x,
                         const double cutoff,
                         const Rcpp::IntegerVector& idx,
                         const int maxit,
                         const double tol,
                         const Rcpp::Nullable<double> th) {
  // parameter dimension
  const int p = par0.size();
  // number of confidence intervals
  const int n = idx.size();
  // confidence intervals (vector)
  std::vector<double> ci_vec;
  ci_vec.reserve(2 * n);
  // #pragma omp parallel for
  for (int j : idx) {
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(j - 1) = 1.0;

    // lower endpoint
    double lower_ub = par0[j - 1];
    double lower_lb = par0[j - 1] - 1.0 / std::log(x.rows());
    // lower bound for lower endpoint
    while (2.0 * MINEL(method, par0, x, lhs,
                       Eigen::Matrix<double, 1, 1>(lower_lb),maxit, tol,
                       th_nloglr(1, th)).nllr <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= 1.0 / std::log(x.rows());
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > tol) {
      if (2.0 * MINEL(method, par0, x, lhs,
                      Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2.0),
                      maxit, tol, th_nloglr(1, th)).nllr > cutoff) {
        lower_lb = (lower_lb + lower_ub) / 2.0;
      } else {
        lower_ub = (lower_lb + lower_ub) / 2.0;
      }
    }
    ci_vec.push_back(lower_ub);

    // upper endpoint
    double upper_lb = par0[j - 1];
    double upper_ub = par0[j - 1] + 1.0 / std::log(x.rows());
    // upper bound for upper endpoint
    while (2.0 * MINEL(method, par0, x, lhs,
                       Eigen::Matrix<double, 1, 1>(upper_ub), maxit, tol,
                       th_nloglr(1, th)).nllr <= cutoff) {
      upper_lb = upper_ub;
      upper_ub += 1.0 / std::log(x.rows());
    }
    // approximate upper bound by numerical search
    while (upper_ub - upper_lb > tol) {
      if (2.0 * MINEL(method, par0, x, lhs,
                      Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2.0),
                      maxit, tol, th_nloglr(1, th)).nllr > cutoff) {
        upper_ub = (upper_lb + upper_ub) / 2.0;
      } else {
        upper_lb = (upper_lb + upper_ub) / 2.0;
      }
    }
    ci_vec.push_back(upper_lb);
  }

  // transform into a matrix
  Eigen::MatrixXd ci = Eigen::Map<Eigen::MatrixXd>(ci_vec.data(), 2, n);
  return ci.transpose();
}

// [[Rcpp::export]]
Eigen::MatrixXd confint_w_(const std::string method,
                           const Eigen::Map<Eigen::VectorXd>& par0,
                           const Eigen::Map<Eigen::MatrixXd>& x,
                           const Eigen::Map<Eigen::ArrayXd>& w,
                           const double cutoff,
                           const Rcpp::IntegerVector& idx,
                           const int maxit,
                           const double tol,
                           const Rcpp::Nullable<double> th) {
  // parameter dimension
  const int p = par0.size();
  // number of confidence intervals
  const int n = idx.size();
  // confidence intervals (vector)
  std::vector<double> ci_vec;
  ci_vec.reserve(2 * n);
  for (int j : idx) {
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(j - 1) = 1.0;

    // lower endpoint
    double lower_ub = par0[j - 1];
    double lower_lb = par0[j - 1] - 1.0 / std::log(x.rows());
    // lower bound for lower endpoint
    while (2.0 * MINEL(method, par0, x, w, lhs,
                       Eigen::Matrix<double, 1, 1>(lower_lb), maxit, tol,
                       th_nloglr(1, th)).nllr <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= 1.0 / std::log(x.rows());
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > tol) {
      if (2.0 * MINEL(method, par0, x, w, lhs,
                      Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2.0),
                      maxit, tol, th_nloglr(1, th)).nllr > cutoff) {
        lower_lb = (lower_lb + lower_ub) / 2.0;
      } else {
        lower_ub = (lower_lb + lower_ub) / 2.0;
      }
    }
    ci_vec.push_back(lower_ub);

    // upper endpoint
    double upper_lb = par0[j - 1];
    double upper_ub = par0[j - 1] + 1.0 / std::log(x.rows());
    // upper bound for upper endpoint
    while (2.0 * MINEL(method, par0, x, w, lhs,
                       Eigen::Matrix<double, 1, 1>(upper_ub), maxit, tol,
                       th_nloglr(1, th)).nllr <= cutoff) {
      upper_lb = upper_ub;
      upper_ub += 1.0 / std::log(x.rows());
    }
    // approximate upper bound by numerical search
    while (upper_ub - upper_lb > tol) {
      if (2.0 * MINEL(method, par0, x, w, lhs,
                      Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2.0),
                      maxit, tol, th_nloglr(1, th)).nllr > cutoff) {
        upper_ub = (upper_lb + upper_ub) / 2.0;
      } else {
        upper_lb = (upper_lb + upper_ub) / 2.0;
      }
    }
    ci_vec.push_back(upper_lb);
  }

  // transform into a matrix
  Eigen::MatrixXd ci = Eigen::Map<Eigen::MatrixXd>(ci_vec.data(), 2, n);
  return ci.transpose();
}
