#include "EL.h"

// [[Rcpp::export]]
Rcpp::List lht_(const std::string method,
                const Eigen::Map<Eigen::VectorXd>& par0,
                const Eigen::Map<Eigen::MatrixXd>& x,
                const Eigen::Map<Eigen::MatrixXd>& lhs,
                const Eigen::Map<Eigen::VectorXd>& rhs,
                const int maxit,
                const double tol,
                const Rcpp::Nullable<double> th) {
  const int q = lhs.rows();
  const int p = lhs.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(lhs);
  if (lu_decomp.rank() != q) {
    Rcpp::stop("'lhs' must have full row rank");
  }

  const MINEL el(method, par0, x, lhs, rhs, maxit, tol, th_nloglr(q, th));
  const double chisq_val = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_val, Rcpp::Named("df") = q,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("method") = method,
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("npar") = p,
    Rcpp::Named("log.prob") = el.logp(x),
    Rcpp::Named("loglik") = el.loglik(),
    Rcpp::Named("coefficients") = el.par,
    Rcpp::Named("statistic") = chisq_val,
    Rcpp::Named("df") = q,
    Rcpp::Named("p.value") = pval);
  return result;
}

// [[Rcpp::export]]
Rcpp::List lht_w_(const std::string method,
                  const Eigen::Map<Eigen::VectorXd>& par0,
                  const Eigen::Map<Eigen::MatrixXd>& x,
                  const Eigen::Map<Eigen::ArrayXd>& w,
                  const Eigen::Map<Eigen::MatrixXd>& lhs,
                  const Eigen::Map<Eigen::VectorXd>& rhs,
                  const int maxit,
                  const double tol,
                  const Rcpp::Nullable<double> th) {
  const int q = lhs.rows();
  const int p = lhs.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(lhs);
  if (lu_decomp.rank() != q) {
    Rcpp::stop("'lhs' must have full row rank");
  }

  const MINEL el(method, par0, x, w, lhs, rhs, maxit, tol, th_nloglr(q, th));
  const double chisq_val = 2.0 * el.nllr;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_val, Rcpp::Named("df") = q,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("method") = method,
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv),
    Rcpp::Named("npar") = p,
    Rcpp::Named("log.prob") = el.logp(x, w),
    Rcpp::Named("loglik") = el.loglik(w),
    Rcpp::Named("coefficients") = el.par,
    Rcpp::Named("statistic") = chisq_val,
    Rcpp::Named("df") = q,
    Rcpp::Named("p.value") = pval);
  return result;
}
