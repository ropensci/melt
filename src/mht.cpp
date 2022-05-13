#include "mht-utils.h"

// [[Rcpp::export]]
Rcpp::List mht_(
    const std::string method,
    const Eigen::Map<Eigen::MatrixXd>& lhs,
    const Eigen::Map<Eigen::MatrixXd>& x,
    const Eigen::Map<Eigen::VectorXd>& est,
    const Eigen::Map<Eigen::VectorXi>& q,
    const int m,
    const int B)
{
  const int p = est.size();
  const Eigen::MatrixXd w = dg0_inv(method, x);
  const Eigen::MatrixXd s = cov(method, est, x);
  // get the square root matrix of the covariance matrix
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
  const Eigen::MatrixXd sqrt_s = es.operatorSqrt();

  Eigen::MatrixXd amat(p, p * m);
  for (int j = 0; j < m; ++j) {
    const Eigen::MatrixXd l = lhs.middleRows(q(j), q(j + 1) - q(j));
    amat.middleCols(j * p, p) = ahat(l, w, s);
  }

  Rcpp::NumericVector out(B);
  for (int b = 0; b < B; ++b) {
    Eigen::RowVectorXd u = rmvn(sqrt_s);
    Eigen::MatrixXd tmp1 = u * amat;
    tmp1.resize(p, m);
    out(b) = (u * tmp1).maxCoeff();
  }

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("s") = s,
    Rcpp::Named("sqrt_s") = sqrt_s,
    Rcpp::Named("amat") = amat,
    Rcpp::Named("out") = out,
    Rcpp::Named("quantile") = quantileRcpp(out, 0.95));
  return result;
}
