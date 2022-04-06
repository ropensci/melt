#include "EL.h"

// [[Rcpp::export]]
Rcpp::List mean_(
    const Eigen::Map<Eigen::VectorXd>& par,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const int maxit,
    const double tol,
    const Rcpp::Nullable<double> th,
    const Rcpp::Nullable<const Eigen::Map<Eigen::ArrayXd>&> w = R_NilValue)
{
  const int n = x.rows();
  const int p = x.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'x' must have full column rank");
  }

  Eigen::VectorXd estimate(p);
  Eigen::ArrayXd w2;
  if (w.isNotNull()) {
    w2 = Rcpp::as<Eigen::ArrayXd>(w);
    estimate = (w2.matrix().transpose() * x) / n;
  } else {
    estimate = x.colwise().mean();
  }

  const EL el = (w.isNull())?
    EL("mean", par, x, maxit, tol, th_nloglr(p, th)) :
    EL("mean", par, x, w2, maxit, tol, th_nloglr(p, th));
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
    Rcpp::Named("log.prob") = el.logp2(x, w),
    Rcpp::Named("loglik") = el.loglik2(w),
    Rcpp::Named("coefficients") = estimate,
    Rcpp::Named("statistic") = chisq_val,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval);
  return result;
}
