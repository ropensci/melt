#include "EL.h"

// [[Rcpp::export]]
Rcpp::NumericVector confreg_(
    const std::string method,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const double cutoff,
    const Eigen::Map<Eigen::MatrixXd>& d,
    const int maxit,
    const int maxit_l,
    const double tol,
    const double tol_l,
    const Rcpp::Nullable<double> step,
    const Rcpp::Nullable<double> th,
    const int nthreads,
    const Eigen::Map<Eigen::ArrayXd>& w)
{
  // confidence region (vector)
  Rcpp::NumericVector reg(d.cols());
  // Eigen::VectorXd reg(d.cols());
  // step size
  const double gamma = step_nloglr(x.rows(), step);
  // test threshold
  const double test_th = th_nloglr(1, th);
  #pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < d.cols(); ++i) {
    Eigen::Matrix2d lhs = Eigen::Matrix2d::Identity();
    const Eigen::Vector2d direction = d.col(i);

    // lower endpoint
    double lower_lb = -1.0;
    double lower_ub = 0;
    // lower bound for lower endpoint
    while (2.0 * MINEL(method, par0, x, lhs, par0 + lower_lb * direction, maxit,
                       maxit_l, tol, tol_l, gamma, test_th, w).nllr <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= 1.0;
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > tol) {
      const double avg = (lower_lb + lower_ub) / 2.0;
      if (2.0 * MINEL(method, par0, x, lhs, par0 + avg * direction, maxit,
                      maxit_l, tol, tol_l, gamma, test_th, w).nllr > cutoff) {
        lower_lb = avg;
      } else {
        lower_ub = avg;
      }
    }
    reg(i) = lower_ub;

    // // upper endpoint
    // double upper_lb = 0;
    // double upper_ub = 1.0;
    // // upper bound for upper endpoint
    // while (2.0 * MINEL(method, par0, x, lhs, par0 + upper_ub * direction, maxit,
    //                    maxit_l, tol, tol_l, gamma, test_th, w).nllr <= cutoff) {
    //   upper_lb = upper_ub;
    //   upper_ub += 1.0;
    // }
    // // approximate upper bound by numerical search
    // while (upper_ub - upper_lb > tol) {
    //   const double avg = (upper_lb + upper_ub) / 2.0;
    //   if (2.0 * MINEL(method, par0, x, lhs, par0 + avg * direction, maxit,
    //                   maxit_l, tol, tol_l, gamma, test_th, w).nllr > cutoff) {
    //     upper_ub = avg;
    //   } else {
    //     upper_lb = avg;
    //   }
    // }
    // reg(i, 1) = upper_lb;
  }
  return reg;
}
