#include "EL.h"

// [[Rcpp::export]]
Rcpp::List lm_(
    const Eigen::Map<Eigen::MatrixXd>& data,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const bool intercept,
    const int maxit,
    const int maxit_l,
    const double tol,
    const double tol_l,
    const Rcpp::Nullable<double> th,
    const int nthreads,
    const Rcpp::Nullable<const Eigen::Map<const Eigen::ArrayXd>&> wt =
      R_NilValue)
{
  const Eigen::VectorXd y = data.col(0);
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const int p = x.cols();

  // overall test
  Eigen::VectorXd par(p);
  Eigen::VectorXd l(p);
  double nllr{};
  int iter{};
  bool conv{};
  Eigen::ArrayXd logp(x.rows());
  double loglik{};
  if (intercept && p > 1) {
    Eigen::MatrixXd lhs(p - 1, p);
    lhs.col(0) = Eigen::MatrixXd::Zero(p - 1, 1);
    lhs.rightCols(p - 1) = Eigen::MatrixXd::Identity(p - 1, p - 1);
    const Eigen::VectorXd rhs = Eigen::VectorXd::Zero(p - 1);
    const MINEL el("lm", par0, data, lhs, rhs, maxit, tol, th_nloglr(p - 1, th),
                   wt);
    par = el.par;
    l = el.l;
    nllr = el.nllr;
    iter = el.iter;
    conv = el.conv;
    logp = el.logp(data);
    loglik = el.loglik();
  } else {
    par = Eigen::VectorXd::Zero(p);
    const EL el("lm", par0, data, maxit_l, tol_l, th_nloglr(p, th), wt);
    l = el.l;
    nllr = el.nllr;
    iter = el.iter;
    conv = el.conv;
    logp = el.logp(data);
    loglik = el.loglik();
  }

  // parameter tests
  Rcpp::NumericVector chisq_val(p);
  Rcpp::LogicalVector par_conv(p);
  Rcpp::Function pchisq("pchisq");
  Rcpp::NumericVector pval(p);
  for (int i = 0; i < p; ++i) {
    Rcpp::checkUserInterrupt();
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(i) = 1.0;
    const MINEL par_test("lm", par0, data, lhs, Eigen::VectorXd::Zero(1),
                         maxit, tol, th_nloglr(1, th), wt);
    chisq_val[i] = 2.0 * par_test.nllr;
    par_conv[i] = par_test.conv;
    pval[i] = Rcpp::as<double>(pchisq(chisq_val[i], Rcpp::Named("df") = 1,
                                      Rcpp::Named("lower.tail") = false));
  }

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("method") = "lm",
      Rcpp::Named("par") = par,
      Rcpp::Named("lambda") = l,
      Rcpp::Named("logLR") = -nllr,
      Rcpp::Named("iterations") = iter,
      Rcpp::Named("convergence") = conv,
      Rcpp::Named("par.tests") = Rcpp::List::create(
        Rcpp::Named("statistic") = chisq_val,
        Rcpp::Named("p.value") = pval,
        Rcpp::Named("convergence") = par_conv)),
    Rcpp::Named("log.prob") = logp,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("statistic") = 2.0 * nllr);
  return result;
}
