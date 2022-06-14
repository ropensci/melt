#include "EL.h"
#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
Rcpp::NumericVector confreg_(
    const std::string method,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const int npar,
    const double cutoff,
    const Rcpp::IntegerVector& idx,
    const Eigen::Map<Eigen::MatrixXd>& circ,
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
  Rcpp::NumericVector cr(circ.cols());
  // step size
  const double gamma = step_nloglr(x.rows(), step);
  // test threshold
  const double test_th = th_nloglr(1, th);
  // lhs matrix
  Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(2, npar);
  lhs(0, idx[0] - 1) = 1;
  lhs(1, idx[1] - 1) = 1;
  // estimates
  const Eigen::Vector2d est{par0(idx[0] - 1), par0(idx[1] - 1)};
  #pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < circ.cols(); ++i) {
    const Eigen::Vector2d direction = circ.col(i);
    // lower endpoint
    double lower_lb = -1.0;
    double lower_ub = 0;
    // lower bound for lower endpoint
    while (2.0 * MINEL(method, par0, x, lhs, est + lower_lb * direction, maxit,
                       maxit_l, tol, tol_l, gamma, test_th, w).nllr <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= 1.0;
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > tol) {
      const double avg = (lower_lb + lower_ub) / 2.0;
      if (2.0 * MINEL(method, par0, x, lhs, est + avg * direction, maxit,
                      maxit_l, tol, tol_l, gamma, test_th, w).nllr > cutoff) {
        lower_lb = avg;
      } else {
        lower_ub = avg;
      }
    }
    cr(i) = lower_ub;
  }
  return cr;
}
