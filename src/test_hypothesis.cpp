#include "utils.h"
#include "EL.h"
#include <RcppEigen.h>
#include <string>

// [[Rcpp::export]]
Rcpp::List test_hypothesis(const std::string method,
                           const Eigen::Map<Eigen::VectorXd> &par0,
                           const Eigen::Map<Eigen::MatrixXd> &x,
                           const Eigen::Map<Eigen::MatrixXd> &lhs,
                           const Eigen::Map<Eigen::VectorXd> &rhs,
                           const int maxit,
                           const int maxit_l,
                           const double tol,
                           const double tol_l,
                           const Rcpp::Nullable<double> step,
                           const Rcpp::Nullable<double> th,
                           const Eigen::Map<Eigen::ArrayXd> &w)
{
    const double test_th = set_threshold(lhs.rows(), th);
    const double gamma = set_step(x.rows(), step);
    const CEL el(method, par0, x, lhs, rhs, maxit, maxit_l, tol, tol_l, gamma,
                 test_th, w);

    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("optim") = Rcpp::List::create(
            Rcpp::Named("par") = el.par,
            Rcpp::Named("lambda") = el.l,
            Rcpp::Named("iterations") = el.iter,
            Rcpp::Named("convergence") = el.conv),
        Rcpp::Named("logp") = el.logp(x, w),
        Rcpp::Named("logl") = el.loglik(w),
        Rcpp::Named("loglr") = -el.nllr,
        Rcpp::Named("statistic") = 2.0 * el.nllr);
    return result;
}
