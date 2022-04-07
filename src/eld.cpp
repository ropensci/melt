#include "EL.h"

// [[Rcpp::export]]
Rcpp::NumericVector eld_(
    const std::string method,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const int maxit,
    const double tol,
    const Rcpp::Nullable<double> th,
    const Rcpp::Nullable<const Eigen::Map<const Eigen::ArrayXd>&> wt =
      R_NilValue)
{
  const int n = x.rows();
  const int p = par0.size();
  const EL el(method, par0, x, maxit, tol, th_nloglr(p, th), wt);
  // maximum empirical likelihood value
  const double mel = el.loglik();
  const Eigen::ArrayXd w = el.w;
  const bool weighted = w.size() != 0;

  Rcpp::NumericVector eld(n);
  // #pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    // Rcpp::checkUserInterrupt();
    // leave-one-out matrix
    Eigen::MatrixXd loo_x(n - 1, x.cols());
    loo_x << x.topRows(i), x.bottomRows(n - 1 - i);
    // leave-one-out weight array
    Eigen::ArrayXd loo_w{};
    if (weighted) {
      loo_w.resize(n - 1);
      loo_w << w.topRows(i), w.bottomRows(n - 1 - i);
      loo_w = ((n - 1) / loo_w.sum()) * loo_w;
    }
    const EL loo_el(method, el.mele_fn(loo_x, loo_w), x, maxit, tol,
                    th_nloglr(p, th), wt);
    eld[i] = 2 * n * (mel - loo_el.loglik());
  }

  return eld;
}
