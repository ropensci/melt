#include "utils.h"

double th_nloglr(const int p, const Rcpp::Nullable<double> th)
{
  return (th.isNull())? 200.0 * p : Rcpp::as<double>(th);
}

Eigen::VectorXd mele_mean(const Eigen::Ref<const Eigen::MatrixXd>& x,
                          const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  if (w.size() == 0) {
    return x.colwise().mean();
  } else {
    return (w.matrix().transpose() * x) / x.rows();
  }
};

Eigen::VectorXd mele_lm(const Eigen::Ref<const Eigen::MatrixXd>& data,
                        const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  const Eigen::VectorXd y = data.col(0);
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  if (w.size() == 0) {
    return x.colPivHouseholderQr().solve(y);
  } else {
    const Eigen::MatrixXd wsqrt =
      Eigen::MatrixXd(w.sqrt().matrix().asDiagonal());
    return (wsqrt * x).colPivHouseholderQr().solve(wsqrt * y);
  }
};


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
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  const int n = g.rows();
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  if (w.size() == 0) {
    return -(1.0 / denominator).sum() * l / n;
  } else {
    return -(w / denominator).sum() * l / n;
  }
}

Eigen::VectorXd gr_nloglr_lm(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  const int n = g.rows();
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  if (w.size() == 0) {
    return
    -(x.transpose() * (x.array().colwise() / denominator).matrix()) * l / n;
  } else {
    return -(x.transpose() * (x.array().colwise() *
             (w / denominator)).matrix()) * l / n;
  }
}







Eigen::ArrayXd logit_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x)
{
  return 1.0 / (1.0 + exp(-x.array()));
}
Eigen::MatrixXd g_bin_logit(const Eigen::Ref<const Eigen::MatrixXd>& data,
                            const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const Eigen::ArrayXd y = data.col(0);
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  return x.array().colwise() * (y - logit_linkinv(x * par));
  // return data.rightCols(data.cols() - 1).array().colwise() *
  //   (data.array().col(0) -
  //   logit_linkinv(data.rightCols(data.cols() - 1) * par));
}
Eigen::VectorXd gr_nloglr_bin_logit(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  const int n = g.rows();
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd numerator =
    logit_linkinv(x * par) * (1.0 - logit_linkinv(x * par));
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  return -(x.transpose() *
           (x.array().colwise() * numerator / denominator).matrix()) * l / n;
}



Eigen::ArrayXd probit_linkinv(const Eigen::Ref<const Eigen::VectorXd>& x)
{
  Eigen::ArrayXd out(x.size());
  for (unsigned int i = 0; i < x.size(); ++i) {
    out[i] = 0.5 * std::erfc(-x[i] * M_SQRT1_2);
  }
  return out;
}
Eigen::MatrixXd g_bin_probit(const Eigen::Ref<const Eigen::MatrixXd>& data,
                             const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const Eigen::ArrayXd y = data.col(0);
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  return x.array().colwise() * (y - probit_linkinv(x * par));
}
Eigen::VectorXd gr_nloglr_bin_probit(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  const int n = g.rows();
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd numerator =
    -exp(-(x * par).array().square() * 0.5) * M_SQRT1_2 * M_2_SQRTPI * 0.5;
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  return -(x.transpose() *
           (x.array().colwise() * numerator / denominator).matrix()) * l / n;
}

















Eigen::ArrayXd linkinv_log(const Eigen::Ref<const Eigen::VectorXd>& x)
{
  return exp(x.array());
}
Eigen::MatrixXd g_gaussian_log(const Eigen::Ref<const Eigen::MatrixXd>& data,
                               const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const Eigen::ArrayXd y = data.col(0);
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  return x.array().colwise() * (y - linkinv_log(x * par));
}
Eigen::VectorXd gr_nloglr_gaussian_log(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const int n = g.rows();
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd numerator = linkinv_log(x * par);
  const Eigen::ArrayXd denominator = Eigen::VectorXd::Ones(n) + g * l;
  return -(x.transpose() *
           (x.array().colwise() * numerator / denominator).matrix()) * l / n;
}

