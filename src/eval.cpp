#include "EL.h"

// [[Rcpp::export]]
Rcpp::List EL_eval(const Eigen::Map<Eigen::MatrixXd>& g,
                   const int maxit,
                   const double abstol,
                   const Rcpp::Nullable<double> threshold)
{
  const int n = g.rows();
  const int p = g.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(g);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'g' must have full column rank");
  }

  const EL el(g, maxit, abstol, th_nlogLR(p, threshold));
  const double chisq_statistic = 2.0 * el.nlogLR;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_statistic, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("lambda") = el.lambda,
      Rcpp::Named("logLR") = -el.nlogLR,
      Rcpp::Named("iterations") = el.iterations,
      Rcpp::Named("convergence") = el.convergence),
        Rcpp::Named("statistic") = chisq_statistic,
        Rcpp::Named("df") = p,
        Rcpp::Named("p.value") = pval,
        Rcpp::Named("alternative") = "two.sided",
        Rcpp::Named("method") = "One sample EL test");
  return result;
}

// [[Rcpp::export]]
Rcpp::List WEL_eval(const Eigen::Map<Eigen::MatrixXd>& g,
                    const Eigen::Map<Eigen::ArrayXd>& w,
                    const int maxit,
                    const double abstol,
                    const Rcpp::Nullable<double> threshold)
{
  const int n = g.rows();
  const int p = g.cols();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(g);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("'g' must have full column rank");
  }

  const EL el(g, w, maxit, abstol, th_nlogLR(p, threshold));
  const double chisq_statistic = 2.0 * el.nlogLR;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_statistic, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("lambda") = el.lambda,
      // Rcpp::Named("log.prob") =  w.log() -
      //   PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g * el.lambda),
      Rcpp::Named("weights") = w,
      Rcpp::Named("logLR") = -el.nlogLR,
      // Rcpp::Named("logWLR") = -el.nlogLR,
      Rcpp::Named("iterations") = el.iterations,
      Rcpp::Named("convergence") = el.convergence),
        Rcpp::Named("statistic") = chisq_statistic,
        Rcpp::Named("df") = p,
        Rcpp::Named("p.value") = pval);
  return result;
}
