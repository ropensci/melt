#include "EL.h"

// [[Rcpp::export]]
Rcpp::List lht_(
    const std::string method,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const Eigen::Map<Eigen::MatrixXd>& lhs,
    const Eigen::Map<Eigen::VectorXd>& rhs,
    const int maxit,
    const int maxit_l,
    const double tol,
    const double tol_l,
    const Rcpp::Nullable<double> th,
    const Eigen::Map<Eigen::ArrayXd>& wt)
{
  const double test_th = th_nloglr(lhs.rows(), th);
  const MINEL el(method, par0, x, lhs, rhs, maxit, maxit_l, tol, tol_l, test_th,
                 wt);

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("method") = method,
      Rcpp::Named("par") = el.par,
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
      Rcpp::Named("log.prob") = el.logp(x),
      Rcpp::Named("loglik") = el.loglik(),
      Rcpp::Named("coefficients") = el.par,
      Rcpp::Named("statistic") = 2.0 * el.nllr);
  return result;
}
