#include "EL.h"

// [[Rcpp::export]]
Rcpp::List EL_mean(const Eigen::Map<Eigen::VectorXd>& par,
                   const Eigen::Map<Eigen::MatrixXd>& x,
                   const int maxit,
                   const double abstol,
                   const Rcpp::Nullable<double> threshold) {
  // check 'par' and 'x'
  const int n = x.rows();
  const int p = x.cols();
  if (par.size() != p) {
    Rcpp::stop("dimensions of 'par' and 'x' do not match");
  }
  if (n < 2) {
    Rcpp::stop("not enough 'x' observations");
  }
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("matrix 'x' must have full column rank");
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
      Rcpp::Named("convergence") = el.convergence,
      Rcpp::Named("control") = Rcpp::List::create(
        Rcpp::Named("maxit") = maxit,
        Rcpp::Named("abstol") = abstol,
        Rcpp::Named("threshold") = th_nlogLR(p, threshold))),
    Rcpp::Named("statistic") = 2.0 * el.nlogLR,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval,
    Rcpp::Named("coefficients") = estimate,
    Rcpp::Named("null.value") = par,
    Rcpp::Named("alternative") = "two.sided",
    Rcpp::Named("method") = "One sample EL test");
  result.attr("class") = Rcpp::CharacterVector({"el_test"});
  return result;
}

