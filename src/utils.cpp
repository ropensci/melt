#include "utils.h"

double th_nloglr(const int p, const Rcpp::Nullable<double> th)
{
  return (th.isNull())? 20.0 * p : Rcpp::as<double>(th);
}

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::MatrixXd>& x,
                       const Eigen::Ref<const Eigen::VectorXd>& par)
{
  return x.rowwise() - par.transpose();
}

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::MatrixXd>& x,
                     const Eigen::Ref<const Eigen::VectorXd>& par)
{
  // const Eigen::VectorXd y = data.col(0);
  // const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  // return x.array().colwise() * (y - x * beta).array();
  return x.rightCols(x.cols() - 1).array().colwise() *
    (x.col(0) - x.rightCols(x.cols() - 1) * par).array();
}

Eigen::VectorXd gr_nloglr_lm(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data)
{
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd denominator =
    Eigen::VectorXd::Ones(g.rows()) + g * lambda;
  return
    -(x.transpose() * (x.array().colwise() / denominator).matrix()) * lambda;
}


