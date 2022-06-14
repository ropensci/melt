#include "EL.h"
#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
Rcpp::List glm_(
    const std::string method,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const Eigen::Map<Eigen::VectorXd>& par0,
    const bool intercept,
    const int maxit,
    const int maxit_l,
    const double tol,
    const double tol_l,
    const Rcpp::Nullable<double> step,
    const Rcpp::Nullable<double> th,
    const int nthreads,
    const Eigen::Map<Eigen::ArrayXd>& w)
{
  const int p = x.cols() - 1;
  const double gamma = step_nloglr(x.rows(), step);

  // overall test
  Eigen::VectorXd par(p);
  Eigen::VectorXd l(p);
  double nllr{};
  int iter{};
  bool conv{};
  Eigen::ArrayXd logp(x.rows());
  double logl{};
  if (intercept && p > 1) {
    Eigen::MatrixXd lhs(p - 1, p);
    lhs.col(0) = Eigen::MatrixXd::Zero(p - 1, 1);
    lhs.rightCols(p - 1) = Eigen::MatrixXd::Identity(p - 1, p - 1);
    const Eigen::VectorXd rhs = Eigen::VectorXd::Zero(p - 1);
    const double test_th = th_nloglr(p - 1, th);
    const MINEL el(method, par0, x, lhs, rhs, maxit, maxit_l, tol, tol_l, gamma,
                   test_th, w);
    par = el.par;
    l = el.l;
    nllr = el.nllr;
    iter = el.iter;
    conv = el.conv;
    logp = el.logp(x, w);
    logl = el.loglik(w);
  } else {
    par = Eigen::VectorXd::Zero(p);
    const double test_th = th_nloglr(p, th);
    const EL el(method, par, x, maxit_l, tol_l, test_th, w);
    l = el.l;
    nllr = el.nllr;
    iter = el.iter;
    conv = el.conv;
    logp = el.logp(x);
    logl = el.loglik();
  }

  // parameter tests
  Rcpp::NumericVector chisq_val(p);
  Rcpp::LogicalVector par_conv(p);
  const double test_th = th_nloglr(1, th);
  // default(none) shared(p, maxit) schedule(auto)
  #pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < p; ++i) {
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(i) = 1.0;
    const MINEL par_test(method, par0, x, lhs, Eigen::VectorXd::Zero(1), maxit,
                         maxit_l, tol, tol_l, gamma, test_th, w);
    chisq_val[i] = 2.0 * par_test.nllr;
    par_conv[i] = par_test.conv;
  }

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("parTests") = Rcpp::List::create(
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










// // [[Rcpp::export]]
// Rcpp::List glm2_(
//     const std::string family,
//     const std::string link,
//     const Eigen::Map<Eigen::MatrixXd>& x,
//     const Eigen::Map<Eigen::VectorXd>& par0,
//     const bool intercept,
//     const int maxit,
//     const int maxit_l,
//     const double tol,
//     const double tol_l,
//     const Rcpp::Nullable<double> step,
//     const Rcpp::Nullable<double> th,
//     const int nthreads,
//     const Eigen::Map<Eigen::ArrayXd>& w)
// {
//   const int n = x.rows();
//   const std::string method = family + "_" + link;
//   const int p = x.cols() - 1;
//   const double gamma = step_nloglr(x.rows(), step);
//
//   Eigen::VectorXd par = par0;
//   par(1) = 0;
//
//   const double test_th = th_nloglr(p - 1, th);
//   const EL el(method, par, x, maxit_l, tol_l, test_th, w);
//
//   const Eigen::MatrixXd g = g_qbin_logit(x, par);
//   Eigen::VectorXd l = (g.transpose() * g).ldlt().solve(g.colwise().sum());
//   double th2 = th_nloglr(p, th);
//   double nllr{0};             // negative log-likelihood ratio
//   int iter{0};                // iterations performed in optimization
//   bool conv{false};           // convergence status
//
//   const PSEUDO_LOG pl(Eigen::VectorXd::Ones(n) + g * l, w);
//   // J matrix
//   Eigen::MatrixXd J;
//
//
//   while (!conv && iter != maxit_l && nllr <= th2) {
//     // pseudo log
//     const PSEUDO_LOG pl(Eigen::VectorXd::Ones(n) + g * l, w);
//     // J matrix
//     J = g.array().colwise() * pl.sqrt_neg_d2plog;
//     // propose new lambda by NR method with least square
//     // Eigen::VectorXd step = (J.transpose() * J).ldlt().solve(
//     //   J.transpose() * (pl.dplog / pl.sqrt_neg_d2plog).matrix());
//     Eigen::VectorXd step =
//       J.colPivHouseholderQr().solve((pl.dplog / pl.sqrt_neg_d2plog).matrix());
//
//     // update function value
//     nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step), w);
//     Rcpp::Rcout << nllr << "\n";
//     // step halving to ensure increase in function value
//     if (nllr < pl.plog_sum) {
//       step /= 2;
//       nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step), w);
//     }
//     // convergence check
//     if (step.norm() < tol_l * l.norm() + tol_l * tol_l) {
//       conv = true;
//     }
//     ++iter;
//     // update lambda
//     l += step;
//   }
//   Eigen::MatrixXd H = J.transpose() * J;
//
//   // // overall test
//   // Eigen::VectorXd par(p + 1);
//   // Eigen::VectorXd l(p + 1);
//   // double nllr{};
//   // int iter{};
//   // bool conv{};
//   // Eigen::ArrayXd logp(x.rows());
//   // double logl{};
//   // if (intercept && p > 1) {
//   //   Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(p - 1, p + 1);
//   //   lhs.middleCols(1, p - 1) = Eigen::MatrixXd::Identity(p - 1, p - 1);
//   //   const Eigen::VectorXd rhs = Eigen::VectorXd::Zero(p - 1);
//   //   const double test_th = th_nloglr(p - 1, th);
//   //   const MINEL el(method, par0, x, lhs, rhs, maxit, maxit_l, tol, tol_l, gamma,
//   //                  test_th, w);
//   //   par = el.par;
//   //   l = el.l;
//   //   nllr = el.nllr;
//   //   iter = el.iter;
//   //   conv = el.conv;
//   //   logp = el.logp(x, w);
//   //   logl = el.loglik(w);
//   // } else {
//   //   Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(p, p + 1);;
//   //   lhs.leftCols(p) = Eigen::MatrixXd::Identity(p, p);
//   //   const Eigen::VectorXd rhs = Eigen::VectorXd::Zero(p);
//   //   const double test_th = th_nloglr(p, th);
//   //   const MINEL el(method, par0, x, lhs, rhs, maxit, maxit_l, tol, tol_l, gamma,
//   //                  test_th, w);
//   //   par = el.par;
//   //   l = el.l;
//   //   nllr = el.nllr;
//   //   iter = el.iter;
//   //   conv = el.conv;
//   //   logp = el.logp(x, w);
//   //   logl = el.loglik(w);
//   // }
//   //
//   // Rcpp::List result = Rcpp::List::create(
//   //   Rcpp::Named("optim") = Rcpp::List::create(
//   //     Rcpp::Named("par") = par,
//   //     Rcpp::Named("lambda") = l,
//   //     Rcpp::Named("iterations") = iter,
//   //     Rcpp::Named("convergence") = conv),
//   //     // Rcpp::Named("parTests") = Rcpp::List::create(
//   //     //   Rcpp::Named("statistic") = chisq_val,
//   //     //   Rcpp::Named("convergence") = par_conv),
//   //   Rcpp::Named("logp") = logp,
//   //   Rcpp::Named("logl") = logl,
//   //   Rcpp::Named("loglr") = -nllr,
//   //   Rcpp::Named("statistic") = 2.0 * nllr
//   //   );
//
//   Rcpp::List result = Rcpp::List::create(
//     Rcpp::Named("par") = par,
//     Rcpp::Named("lambda") = l,
//     Rcpp::Named("lambda2") = el.l,
//     Rcpp::Named("iterations") = iter,
//     Rcpp::Named("convergence") = conv,
//     Rcpp::Named("nllr") = nllr,
//     Rcpp::Named("nllr2") = el.nllr,
//     Rcpp::Named("tmp") = g_qbin_logit(x, par),
//     Rcpp::Named("h") = H);
//   return result;
// }
