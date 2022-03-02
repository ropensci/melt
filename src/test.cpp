#include "EL.h"

// [[Rcpp::export]]
Rcpp::List EL_test(const std::string method,
                   const Eigen::Map<Eigen::VectorXd>& par,
                   const Eigen::Map<Eigen::MatrixXd>& x,
                   const int maxit,
                   const double abstol,
                   const Rcpp::Nullable<double> threshold) {
  const int p = par.size();
  const EL2 el(par, x, method, maxit, abstol, th_nlogLR(p, threshold));
  const double chisq_statistic = 2.0 * el.nlogLR;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_statistic, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("type") = method,
      Rcpp::Named("lambda") = el.lambda,
      Rcpp::Named("logLR") = -el.nlogLR,
      Rcpp::Named("iterations") = el.iterations,
      Rcpp::Named("convergence") = el.convergence,
      Rcpp::Named("control") = Rcpp::List::create(
        Rcpp::Named("maxit") = maxit,
        Rcpp::Named("abstol") = abstol,
        Rcpp::Named("threshold") = th_nlogLR(p, threshold))),
        Rcpp::Named("statistic") = chisq_statistic,
        Rcpp::Named("df") = p,
        Rcpp::Named("p.value") = pval,
        Rcpp::Named("null.value") = par,
        Rcpp::Named("alternative") = "two.sided",
        Rcpp::Named("method") = "One sample EL test");
  return result;
}
