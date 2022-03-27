#include "EL.h"

// [[Rcpp::export]]
Rcpp::List lm_(const Eigen::MatrixXd& data,
               const bool intercept,
               const int maxit,
               const double tol,
               const Rcpp::Nullable<double> th) {
  // check design matrix
  const Eigen::VectorXd y = data.col(0);
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const int p = x.cols();
  if (x.rows() <= p) {
    Rcpp::stop("design matrix must have full column rank");
  }
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p) {
    Rcpp::stop("design matrix must have full column rank");
  }

  // overall test
  const Eigen::VectorXd bhat = x.colPivHouseholderQr().solve(y);
  const Eigen::VectorXd fit_val = x * bhat;
  const Eigen::VectorXd resid = y - fit_val;
  Eigen::MatrixXd lhs(p - 1, p);
  lhs.col(0) = Eigen::MatrixXd::Zero(p - 1, 1);
  lhs.rightCols(p - 1) = Eigen::MatrixXd::Identity(p - 1, p - 1);
  const Eigen::VectorXd rhs = Eigen::VectorXd::Zero(p - 1);
  const EL el =
    (intercept && p > 1)?
    EL("lm", bhat, data, lhs, rhs, maxit, tol, th_nloglr(p - 1, th)) :
    EL("lm", Eigen::MatrixXd::Zero(1, p), data , maxit, tol, th_nloglr(p, th));

  // test each coefficient
  Rcpp::NumericVector chisq_val(p);
  Rcpp::LogicalVector conv(p);
  Rcpp::Function pchisq("pchisq");
  Rcpp::NumericVector pval(p);
  for (int i = 0; i < p; ++i) {
    Rcpp::checkUserInterrupt();
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
    lhs(i) = 1.0;
    const EL par_test("lm", bhat, data, lhs, Eigen::VectorXd::Zero(1), maxit,
                       tol, th_nloglr(1, th));
    chisq_val[i] = 2.0 * par_test.nllr;
    conv[i] = par_test.conv;
    pval[i] = Rcpp::as<double>(pchisq(chisq_val[i], Rcpp::Named("df") = 1,
                                      Rcpp::Named("lower.tail") = false));
  }

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("method") = "lm",
      Rcpp::Named("lambda") = el.l,
      Rcpp::Named("logLR") = -el.nllr,
      Rcpp::Named("iterations") = el.iter,
      Rcpp::Named("convergence") = el.conv,
      Rcpp::Named("par.tests") = Rcpp::List::create(
        Rcpp::Named("statistic") = chisq_val,
        Rcpp::Named("p.value") = pval,
        Rcpp::Named("convergence") = conv)),
    Rcpp::Named("coefficients") = bhat,
    Rcpp::Named("residuals") = resid,
    Rcpp::Named("rank") = p,
    Rcpp::Named("fitted.values") = fit_val);
  result.attr("class") = Rcpp::CharacterVector({"el_lm", "el_test"});
  return result;
}
