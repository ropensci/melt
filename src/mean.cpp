#include "EL.h"

// [[Rcpp::export]]
Rcpp::List EL_mean(const Eigen::Map<Eigen::VectorXd>& par,
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
  const double chisq_statistic = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval =
    Rcpp::as<double>(pchisq(chisq_statistic, Rcpp::Named("df") = p,
                            Rcpp::Named("lower.tail") = false));
  const Eigen::VectorXd estimate = x.colwise().mean();

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("type") = "mean",
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("statistic") = chisq_statistic,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval,
    Rcpp::Named("coefficients") = estimate,
    Rcpp::Named("null.value") = par);
  result.attr("class") = Rcpp::CharacterVector({"el_test"});
  return result;
}

// [[Rcpp::export]]
Rcpp::List WEL_mean(const Eigen::Map<Eigen::VectorXd>& par,
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

  EL el("mean", par, x, w, maxit, tol, th_nloglr(p, th));
  const double chisq_statistic = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval =
    Rcpp::as<double>(pchisq(chisq_statistic, Rcpp::Named("df") = p,
                            Rcpp::Named("lower.tail") = false));
  const Eigen::VectorXd estimate = x.colwise().mean();

  Eigen::ArrayXd log_prob = el.log_prob(x, w);
  Eigen::ArrayXd log_wprob = el.log_wprob(x, w);

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("type") = "mean",
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("log.prob") = log_prob,
      Rcpp::Named("log.wprob") = log_wprob,
      Rcpp::Named("weights") = w,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("statistic") = chisq_statistic,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval,
    Rcpp::Named("coefficients") = estimate,
    Rcpp::Named("null.value") = par);
  result.attr("class") = Rcpp::CharacterVector({"el_test"});
  return result;
}
