#include "utils.h"
#include "EL.h"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string>
#include <vector>

// [[Rcpp::export]]
Rcpp::List test_LM(const Eigen::Map<Eigen::MatrixXd> &x,
                   const Eigen::Map<Eigen::VectorXd> &par0,
                   const bool intercept,
                   const int maxit,
                   const int maxit_l,
                   const double tol,
                   const double tol_l,
                   const Rcpp::Nullable<double> step,
                   const Rcpp::Nullable<double> th,
                   const int nthreads,
                   const Eigen::Map<Eigen::ArrayXd> &w)
{
  const int p = x.cols() - 1;
  const double gamma = set_step(x.rows(), step);

  // overall test
  Eigen::VectorXd par(p);
  Eigen::VectorXd l(p);
  double nllr{};
  int iter{};
  bool conv{};
  Eigen::ArrayXd logp(x.rows());
  double logl{};
  if (intercept && p > 1)
  {
    Eigen::MatrixXd lhs(p - 1, p);
    lhs.col(0) = Eigen::MatrixXd::Zero(p - 1, 1);
    lhs.rightCols(p - 1) = Eigen::MatrixXd::Identity(p - 1, p - 1);
    const Eigen::VectorXd rhs = Eigen::VectorXd::Zero(p - 1);
    const double test_th = set_threshold(p - 1, th);
    const CEL el("lm", par0, x, lhs, rhs, maxit, maxit_l, tol, tol_l, gamma,
                 test_th, w);
    par = el.par;
    l = el.l;
    nllr = el.nllr;
    iter = el.iter;
    conv = el.conv;
    logp = el.logp(x, w);
    logl = el.loglik(w);
  }
  else
  {
    par = Eigen::VectorXd::Zero(p);
    const double test_th = set_threshold(p, th);
    const EL el("lm", par, x, maxit_l, tol_l, test_th, w);
    l = el.l;
    nllr = el.nllr;
    iter = el.iter;
    conv = el.conv;
    logp = el.logp(x);
    logl = el.loglik();
  }

  // parameter tests
  std::vector<double> chisq_val(p);
  std::vector<bool> par_conv(p);
  const double test_th = set_threshold(1, th);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
  for (int i = 0; i < p; ++i)
  {
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(i) = 1.0;
    const CEL par_test("lm", par0, x, lhs, Eigen::VectorXd::Zero(1), maxit,
                       maxit_l, tol, tol_l, gamma, test_th, w);
    chisq_val[i] = 2.0 * par_test.nllr;
    par_conv[i] = par_test.conv;
  }

  Rcpp::List result = Rcpp::List::create(
      Rcpp::Named("par_tests") = Rcpp::List::create(
          Rcpp::Named("statistic") = chisq_val,
          Rcpp::Named("convergence") = par_conv),
      Rcpp::Named("optim") = Rcpp::List::create(
          Rcpp::Named("par") = par,
          Rcpp::Named("lambda") = l,
          Rcpp::Named("iterations") = iter,
          Rcpp::Named("convergence") = conv),
      Rcpp::Named("logp") = logp,
      Rcpp::Named("logl") = logl,
      Rcpp::Named("loglr") = -nllr,
      Rcpp::Named("statistic") = 2.0 * nllr);
  return result;
}
