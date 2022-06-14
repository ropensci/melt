#include "EL.h"
#include "utils.h"

// [[Rcpp::export]]
Rcpp::List eval_(
    const std::string method,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const int maxit_l,
    const double tol_l,
    const Rcpp::Nullable<double> th,
    const Eigen::Map<Eigen::ArrayXd>& wt)
{
  const double test_th = th_nloglr(par0.size(), th);
  const EL el(method, par0, x, maxit_l, tol_l, test_th, wt);

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("par") = el.par,
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("logp") = el.logp(x),
    Rcpp::Named("logl") = el.loglik(),
    Rcpp::Named("loglr") = -el.nllr,
    Rcpp::Named("statistic") = 2.0 * el.nllr);
  return result;
}

// [[Rcpp::export]]
Rcpp::List eval_g_(
    const Eigen::Map<Eigen::MatrixXd>& g,
    const int maxit_l,
    const double tol_l,
    const Rcpp::Nullable<double> th,
    const Eigen::Map<Eigen::ArrayXd>& wt)
{
  const double test_th = th_nloglr(g.cols(), th);
  const EL el(g, maxit_l, tol_l, test_th, wt);

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("logp") = el.logp_g(g),
    Rcpp::Named("logl") = el.loglik(),
    Rcpp::Named("loglr") = -el.nllr,
    Rcpp::Named("statistic") = 2.0 * el.nllr);
  return result;
}
