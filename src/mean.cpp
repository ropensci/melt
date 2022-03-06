#include "EL.h"

// [[Rcpp::export]]
Rcpp::List EL_mean(const Eigen::Map<Eigen::VectorXd>& par,
                   const Eigen::Map<Eigen::MatrixXd>& x,
                   const int maxit,
                   const double abstol,
                   const Rcpp::Nullable<double> threshold) {
  const int n = x.rows();
  const int p = x.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'x' must have full column rank");
  }

  const EL2 el(par, x, "mean", maxit, abstol, th_nlogLR(p, threshold));
  const double chisq_statistic = 2.0 * el.nlogLR;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_statistic, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));
  const Eigen::VectorXd estimate = x.colwise().mean();

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("type") = "mean",
      Rcpp::Named("lambda") = el.lambda,
      Rcpp::Named("logLR") = -el.nlogLR,
      Rcpp::Named("iterations") = el.iterations,
      Rcpp::Named("convergence") = el.convergence),
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
                    const double abstol,
                    const Rcpp::Nullable<double> threshold) {
  const int n = x.rows();
  const int p = x.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'x' must have full column rank");
  }

  EL2 el(par, x, w, "mean", maxit, abstol, th_nlogLR(p, threshold));
  const double chisq_statistic = 2.0 * el.nlogLR;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_statistic, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));
  const Eigen::VectorXd estimate = x.colwise().mean();

  const Eigen::ArrayXd log_prob = el.log_prob(x, w);
  const Eigen::ArrayXd log_wprob = el.log_wprob(x, w);

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("type") = "mean",
      Rcpp::Named("lambda") = el.lambda,
      // Rcpp::Named("log.prob") = log_prob,
      // Rcpp::Named("log.wprob") = log_wprob,
      Rcpp::Named("weights") = w,
      // Rcpp::Named("logLR") = log_prob.sum(),
      Rcpp::Named("logLR") = -el.nlogLR,
      Rcpp::Named("iterations") = el.iterations,
      Rcpp::Named("convergence") = el.convergence),
        Rcpp::Named("statistic") = chisq_statistic,
        Rcpp::Named("df") = p,
        Rcpp::Named("p.value") = pval,
        Rcpp::Named("coefficients") = estimate,
        Rcpp::Named("null.value") = par);
  result.attr("class") = Rcpp::CharacterVector({"el_test"});
  return result;
}
