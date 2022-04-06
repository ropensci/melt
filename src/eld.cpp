#include "EL.h"

// [[Rcpp::export]]
Rcpp::NumericVector eld_(
    const std::string method,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const int maxit,
    const double tol,
    const Rcpp::Nullable<double> th,
    const Rcpp::Nullable<const Eigen::Map<Eigen::ArrayXd>&> w = R_NilValue)
{
  const int p = par0.size();
  const int n = x.rows();

  const EL el = (w.isNull())?
    EL(method, par0, x, maxit, tol, th_nloglr(p, th)) :
    EL(method, par0, x, Rcpp::as<Eigen::ArrayXd>(w), maxit, tol,
       th_nloglr(p, th));
  // maximum empirical likelihood value
  const double mel = el.loglik2(w);

  Rcpp::NumericVector eld(n);
  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    // leave-one-out matrix
    Eigen::MatrixXd loo(n - 1, x.cols());
    loo << x.topRows(i), x.bottomRows(n - 1 - i);
    const EL loo_el(method, el.mele_fcn(loo), x, maxit, tol, th_nloglr(p, th));
    eld[i] = 2 * n * (mel - loo_el.loglik());
  }

  // Rcpp::List result = Rcpp::List::create(
  //   Rcpp::Named("eld") = eld);
  return eld;
}
