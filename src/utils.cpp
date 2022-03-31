#include "utils.h"

double th_nloglr(const int p, const Rcpp::Nullable<double> th)
{
  return (th.isNull())? 200.0 * p : Rcpp::as<double>(th);
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

Eigen::VectorXd gr_nloglr_mean(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const int n = g.rows();
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  return -(1.0 / denominator).sum() * l / n;
}

Eigen::VectorXd wgr_nloglr_mean(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  const int n = g.rows();
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  return -(w / denominator).sum() * l / n;
}

Eigen::VectorXd gr_nloglr_lm(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const int n = g.rows();
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  return
    -(x.transpose() * (x.array().colwise() / denominator).matrix()) * l / n;
}

Eigen::VectorXd wgr_nloglr_lm(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  const int n = g.rows();
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  return -(x.transpose() * (x.array().colwise() * (w / denominator)).matrix()) *
    l / n;
}








Eigen::ArrayXd logit_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x)
{
  return 1.0 / (1.0 + exp(-x.array()));
}

Eigen::MatrixXd g_logit(const Eigen::Ref<const Eigen::MatrixXd>& data,
                        const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const Eigen::ArrayXd y = data.col(0);
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  return x.array().colwise() * (y - logit_linkinv(x * par));
  // return data.rightCols(data.cols() - 1).array().colwise() *
  //   (data.array().col(0) -
  //   logit_linkinv(data.rightCols(data.cols() - 1) * par));
}

Eigen::VectorXd gr_nloglr_logit(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const int n = g.rows();
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd numerator =
    logit_linkinv(x * par) * (1.0 - logit_linkinv(x * par));
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  return -(x.transpose() *
           (x.array().colwise() * numerator / denominator).matrix()) * l / n;
}
