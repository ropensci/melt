#include "EL.h"

// [[Rcpp::export]]
Rcpp::List mean_(const Eigen::Map<Eigen::VectorXd>& par,
                 const Eigen::Map<Eigen::MatrixXd>& x,
                 const int maxit,
                 const double tol,
                 const Rcpp::Nullable<double> th) {
  const int n = x.rows();
  const int p = x.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'x' must have full column rank");
  }

  const EL el("mean", par, x, maxit, tol, th_nloglr(p, th));
  const double chisq_val = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval =
    Rcpp::as<double>(pchisq(chisq_val, Rcpp::Named("df") = p,
                            Rcpp::Named("lower.tail") = false));
  const Eigen::VectorXd estimate = x.colwise().mean();

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
  result.attr("class") = Rcpp::CharacterVector({"el_test"});
  return result;
}

// [[Rcpp::export]]
Rcpp::List mean_w_(const Eigen::Map<Eigen::VectorXd>& par,
                   const Eigen::Map<Eigen::MatrixXd>& x,
                   const Eigen::Map<Eigen::ArrayXd>& w,
                   const int maxit,
                   const double tol,
                   const Rcpp::Nullable<double> th) {
  const int n = x.rows();
  const int p = x.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'x' must have full column rank");
  }

  const EL el("mean", par, x, w, maxit, tol, th_nloglr(p, th));
  const double chisq_val = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval =
    Rcpp::as<double>(pchisq(chisq_val, Rcpp::Named("df") = p,
                            Rcpp::Named("lower.tail") = false));
  const Eigen::VectorXd estimate = (w.matrix().transpose() * x) / n;

  // const Eigen::ArrayXd log_prob = el.log_prob(x, w);
  // const Eigen::ArrayXd log_wprob = el.log_wprob(x, w);

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("method") = "mean",
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("npar") = p,
    Rcpp::Named("log.prob") = el.logp(x, w),
    Rcpp::Named("loglik") = el.wloglik(w),
    Rcpp::Named("coefficients") = estimate,
    Rcpp::Named("statistic") = chisq_val,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval);
  result.attr("class") = Rcpp::CharacterVector({"el_test"});
  return result;
}
