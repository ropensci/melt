#include "EL.h"

// [[Rcpp::export]]
Rcpp::List EL_lht(const std::string method,
                  const Eigen::Map<Eigen::VectorXd>& par0,
                  const Eigen::Map<Eigen::MatrixXd>& x,
                  const Eigen::Map<Eigen::MatrixXd>& lhs,
                  const Eigen::Map<Eigen::VectorXd>& rhs,
                  const int maxit,
                  const double tol,
                  const Rcpp::Nullable<double> th) {
  const int p = lhs.rows();
  const EL el(method, par0, x, lhs, rhs, maxit, tol, th_nloglr(p, th));
  const double chisq_val = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_val, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("type") = method,
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("statistic") = chisq_val,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval);
  return result;
}
