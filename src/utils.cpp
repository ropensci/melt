#include "utils.h"

double th_nlogLR(const int p, const Rcpp::Nullable<double> threshold) {
  return (threshold.isNull())? 20.0 * p : Rcpp::as<double>(threshold);
}

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::VectorXd>& par,
                       const Eigen::Ref<const Eigen::MatrixXd>& x) {
  return x.rowwise() - par.transpose();
}

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::VectorXd>& par,
                     const Eigen::Ref<const Eigen::MatrixXd>& data) {
  // const Eigen::VectorXd y = data.col(0);
  // const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  // return x.array().colwise() * (y - x * beta).array();
  return data.rightCols(data.cols() - 1).array().colwise() *
    (data.col(0) - data.rightCols(data.cols() - 1) * par).array();
}

Eigen::VectorXd gr_nlogLR_lm(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data) {
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd denominator =
    Eigen::VectorXd::Ones(g.rows()) + g * lambda;
  // const Eigen::MatrixXd gradient =
  //   -(x.transpose() * (x.array().colwise() / denominator).matrix()) * lambda;
  return
    -(x.transpose() * (x.array().colwise() / denominator).matrix()) * lambda;
}


