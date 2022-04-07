#include "EL.h"

// [[Rcpp::export]]
Rcpp::List mean_(
    const Eigen::Map<Eigen::VectorXd>& par,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const int maxit,
    const double tol,
    const Rcpp::Nullable<double> th,
    const Rcpp::Nullable<const Eigen::Map<const Eigen::ArrayXd>&> wt =
      R_NilValue)
{
  const int n = x.rows();
  const int p = par.size();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'x' must have full column rank");
  }

  Eigen::VectorXd estimate(p);
  if (wt.isNull()) {
    estimate = x.colwise().mean();
  } else {
    estimate = (Rcpp::as<Eigen::VectorXd>(wt).transpose() * x) / n;
  }

  const EL el("mean", par, x, maxit, tol, th_nloglr(p, th), wt);

  const double chisq_val = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval =
    Rcpp::as<double>(pchisq(chisq_val, Rcpp::Named("df") = p,
                            Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("method") = "mean",
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("npar") = p,
    Rcpp::Named("log.prob") = el.logp(x),
    Rcpp::Named("loglik") = el.loglik(),
    Rcpp::Named("coefficients") = estimate,
    Rcpp::Named("statistic") = chisq_val,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval);
  return result;
}
