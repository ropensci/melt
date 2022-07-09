#include "test_multiple_hypotheses_utils.h"
#include "helpers.h"
#include <RcppEigen.h>
#include <functional>
#include <map>
#include <string>
#include <vector>

Eigen::RowVectorXd rmvn(const Eigen::Ref<const Eigen::MatrixXd> &sqrt)
{
  Eigen::RowVectorXd u(sqrt.cols());
  for (int i = 0; i < sqrt.cols(); ++i)
  {
    u(i) = R::rnorm(0, 1.0);
  }
  return u * sqrt;
}

Eigen::MatrixXd w_mean(const Eigen::Ref<const Eigen::MatrixXd> &x,
                       const Eigen::Ref<const Eigen::VectorXd> &par)
{
  return Eigen::MatrixXd::Identity(x.cols(), x.cols());
}

Eigen::MatrixXd w_lm(const Eigen::Ref<const Eigen::MatrixXd> &x,
                     const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  return static_cast<double>(x.rows()) * (xmat.transpose() * xmat).inverse();
}

Eigen::MatrixXd w_gauss_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c = y * log_linkinv(xmat * par) -
                           2.0 * log_linkinv(2.0 * xmat * par);
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_gauss_inverse(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c = cube(inverse_linkinv(xmat * par)) *
                           (2.0 * y - 3.0 * inverse_linkinv(xmat * par));
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_bin_logit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  Eigen::ArrayXd c = logit_linkinv(xmat * par) *
                     (1.0 - logit_linkinv(xmat * par));
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_bin_probit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd phi = 0.5 * M_SQRT1_2 * M_2_SQRTPI *
                             exp(-0.5 * (xmat * par).array().square());
  const Eigen::ArrayXd dphi = -0.5 * M_SQRT1_2 * M_2_SQRTPI *
                              exp(-0.5 * (xmat * par).array().square()) * (xmat * par).array();
  const Eigen::ArrayXd c = y *
                               (dphi * inverse(probit_linkinv(xmat * par)) -
                                square(inverse(probit_linkinv(xmat * par)) * phi)) -
                           dphi;
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_bin_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                          const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c = square((1.0 - log_linkinv(xmat * par)).inverse()) *
                           log_linkinv(xmat * par) * (y - 1.0);
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_poi_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                          const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c = log_linkinv(xmat * par);
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_poi_identity(const Eigen::Ref<const Eigen::MatrixXd> &x,
                               const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c = y * square(inverse((xmat * par).array()));
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_poi_sqrt(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c =
      2.0 * (y * square(((xmat * par).array().inverse()))) + 2.0;
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd dg0_inv(const std::string method,
                        const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::VectorXd> &par)
{
  std::map<std::string, std::function<Eigen::MatrixXd(
                            const Eigen::Ref<const Eigen::MatrixXd> &,
                            const Eigen::Ref<const Eigen::VectorXd> &)>>
      w_map{{{"mean", w_mean},
             {"lm", w_lm},
             {"gaussian_identity", w_lm},
             {"gaussian_log", w_gauss_log},
             {"gaussian_inverse", w_gauss_inverse},
             {"binomial_logit", w_bin_logit},
             {"binomial_probit", w_bin_probit},
             {"binomial_log", w_bin_log},
             {"poisson_log", w_poi_log},
             {"poisson_identity", w_poi_identity},
             {"poisson_sqrt", w_poi_sqrt}}};
  return w_map[method](x, par);
}

Eigen::MatrixXd shat(const std::string method,
                     const Eigen::Ref<const Eigen::VectorXd> &par,
                     const Eigen::Ref<const Eigen::MatrixXd> &x)
{
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd> &,
      const Eigen::Ref<const Eigen::VectorXd> &)>
      g_fn = set_g_fn(method);
  return (1.0 / x.rows()) * ((g_fn(x, par)).transpose() * g_fn(x, par));
}

Eigen::MatrixXd ahat(const Eigen::Ref<const Eigen::MatrixXd> &j,
                     const Eigen::Ref<const Eigen::MatrixXd> &w,
                     const Eigen::Ref<const Eigen::MatrixXd> &s)
{
  return (j * w).transpose() * (((j * w) * s * (j * w).transpose()).inverse()) *
         (j * w);
}
