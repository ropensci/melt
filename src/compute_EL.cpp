#include "utils.h"
#include "EL.h"
#include <RcppEigen.h>
#include <string>

// [[Rcpp::export]]
Rcpp::List compute_EL(const std::string method,
                      const Eigen::Map<Eigen::VectorXd> &par0,
                      const Eigen::Map<Eigen::MatrixXd> &x,
                      const int maxit_l,
                      const double tol_l,
                      const Rcpp::Nullable<double> th,
                      const Eigen::Map<Eigen::ArrayXd> &w)
{
    const double test_th = set_threshold(par0.size(), th);
    const EL el(method, par0, x, maxit_l, tol_l, test_th, w);

    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("optim") = Rcpp::List::create(
            Rcpp::Named("par") = el.par,
            Rcpp::Named("lambda") = el.l,
            Rcpp::Named("iterations") = el.iter,
            Rcpp::Named("convergence") = el.conv),
        Rcpp::Named("logp") = el.logp(x),
        Rcpp::Named("logl") = el.loglik(),
        Rcpp::Named("loglr") = -el.nllr,
        Rcpp::Named("statistic") = 2.0 * el.nllr);
    return result;
}
