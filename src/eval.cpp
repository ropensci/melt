#include "EL.h"

// [[Rcpp::export]]
Rcpp::List eval_(const Eigen::Map<Eigen::MatrixXd>& g,
                 const int maxit,
                 const double tol,
                 const Rcpp::Nullable<double> th)
{
  const int n = g.rows();
  const int p = g.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(g);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'g' must have full column rank");
  }

  const EL el(g, maxit, tol, th_nloglr(p, th));
  const double chisq_val = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_val, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("npar") = p,
    Rcpp::Named("log.prob") = el.logp_g(g),
    Rcpp::Named("loglik") = el.loglik(),
    Rcpp::Named("statistic") = chisq_val,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval);
  return result;
}

// [[Rcpp::export]]
Rcpp::List eval_w_(const Eigen::Map<Eigen::MatrixXd>& g,
                   const Eigen::Map<Eigen::ArrayXd>& w,
                   const int maxit,
                   const double tol,
                   const Rcpp::Nullable<double> th)
{
  const int n = g.rows();
  const int p = g.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(g);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'g' must have full column rank");
  }

  const EL el(g, w, maxit, tol, th_nloglr(p, th));
  const double chisq_val = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_val, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("npar") = p,
    Rcpp::Named("log.prob") = el.logp_g(g, w),
    Rcpp::Named("loglik") = el.wloglik(w),
    Rcpp::Named("statistic") = chisq_val,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval);
  return result;
}
