#include "EL.h"

Rcpp::List tmp_(
    const Eigen::Map<Eigen::MatrixXd>& x,
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
  const int p = x.cols() - 1;
  const int n = x.rows();
  Eigen::ArrayXd w;
  if (wt.isNotNull()) {
    w = Rcpp::as<Eigen::ArrayXd>(wt);
  }
  Eigen::MatrixXd lhs(p - 1, p);
  lhs.col(0) = Eigen::MatrixXd::Zero(p - 1, 1);
  lhs.rightCols(p - 1) = Eigen::MatrixXd::Identity(p - 1, p - 1);
  const Eigen::VectorXd rhs = Eigen::VectorXd::Zero(p - 1);
  const double test_th = th_nloglr(p - 1, th);

  // orthogonal projection matrix
  const Eigen::MatrixXd proj =
    Eigen::MatrixXd::Identity(lhs.cols(), lhs.cols()) -
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * lhs;
  // parameter (constraint imposed)
  Eigen::VectorXd par = proj * par0 + lhs.transpose() * (lhs * lhs.transpose()).inverse() * rhs;
  // estimating function
  Eigen::MatrixXd g = g_lm(x, par);
  // lambda
  Eigen::VectorXd l = EL(g, maxit_l, tol_l, test_th, wt).l;
  // function value (-logLR)
  double nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * l, w);
  // function norm
  // const double norm0 = (proj * gr_nloglr_lm(l, g, x, par, w)).norm();
  //
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  Eigen::MatrixXd xx = xmat.array().colwise() / denominator;
  // Eigen::MatrixXd xxx = xmat.transpose() * (xmat.array().colwise() / denominator);
  // return -(xmat.transpose() * xx) * l;
  // return -(x.transpose() * (x.array().colwise() / denominator).matrix()) * l;

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("proj") = proj,
      Rcpp::Named("par") = par,
      Rcpp::Named("l") = l),
      Rcpp::Named("nllr") = nllr,
      Rcpp::Named("grad") = gr_nloglr_lm(l, g, x, par, w),
      Rcpp::Named("t1") = xx,
      Rcpp::Named("t2") = (xmat.array().colwise() / denominator).matrix(),
      Rcpp::Named("t3") = xmat.transpose() * xx,
      Rcpp::Named("t4") = xmat.transpose() *
        (xmat.array().colwise() / denominator).matrix(),
      Rcpp::Named("t5") =
        -(xmat.transpose() * (xmat.array().colwise() / denominator).matrix()) * l,
      Rcpp::Named("t6") =-(xmat.transpose() * xx) * l);
  return result;
}
