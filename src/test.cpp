#include "EL.h"

// [[Rcpp::export]]
Rcpp::List EL_test(const std::string method,
                   const Eigen::Map<Eigen::VectorXd>& par,
                   const Eigen::Map<Eigen::MatrixXd>& x,
                   const int maxit,
                   const double abstol,
                   const Rcpp::Nullable<double> th) {
  const int p = par.size();
  const EL el(method, par, x, maxit, abstol, th_nloglr(p, th));
  const double chisq_statistic = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_statistic, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("type") = method,
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("statistic") = chisq_statistic,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval,
    Rcpp::Named("null.value") = par,
    Rcpp::Named("alternative") = "two.sided",
    Rcpp::Named("method") = "One sample EL test");
  return result;
}
