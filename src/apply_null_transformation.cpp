#include "apply_null_transformation.h"
#include "helpers.h"
#include <RcppEigen.h>
#include <functional>
#include <map>
#include <string>

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &)>
transform_x_fn(const std::string method)
{
  std::map<std::string, std::function<Eigen::MatrixXd(
                            const Eigen::Ref<const Eigen::MatrixXd> &,
                            const Eigen::Ref<const Eigen::VectorXd> &,
                            const Eigen::Ref<const Eigen::VectorXd> &)>>
      x0_map{{{"mean", x0_mean},
              {"sd", x0_sd},
              {"lm", x0_lm},
              {"gaussian_identity", x0_lm},
              {"gaussian_log", x0_gauss_log},
              {"gaussian_inverse", x0_gauss_inverse},
              {"binomial_logit", x0_bin_logit},
              {"binomial_probit", x0_bin_probit},
              {"binomial_log", x0_bin_log},
              {"poisson_log", x0_poi_log},
              {"poisson_identity", x0_poi_identity},
              {"poisson_sqrt", x0_poi_sqrt},
              {"quasipoisson_log", x0_poi_log},
              {"quasipoisson_identity", x0_poi_identity},
              {"quasipoisson_sqrt", x0_poi_sqrt}}};
  return x0_map[method];
}

Eigen::MatrixXd x0_mean(const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::VectorXd> &par,
                        const Eigen::Ref<const Eigen::VectorXd> &est)
{
  return x.rowwise() + (par - est).transpose();
}

Eigen::MatrixXd x0_sd(const Eigen::Ref<const Eigen::MatrixXd> &x,
                      const Eigen::Ref<const Eigen::VectorXd> &par,
                      const Eigen::Ref<const Eigen::VectorXd> &est)
{
  return x.rowwise() +
         (par.array().square() - est.array().square()).matrix().transpose();
}

Eigen::MatrixXd x0_lm(const Eigen::Ref<const Eigen::MatrixXd> &x,
                      const Eigen::Ref<const Eigen::VectorXd> &par,
                      const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::VectorXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out = x;
  out.col(1) = y + (xmat * (par - est));
  return out;
}

Eigen::MatrixXd x0_gauss_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par,
                             const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out = x;
  out.col(1) =
      log_linkinv(xmat * par + s) +
      (y - log_linkinv(xmat * est + s)) * log_linkinv(xmat * (est - par));
  return out;
}

Eigen::MatrixXd x0_gauss_inverse(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                 const Eigen::Ref<const Eigen::VectorXd> &par,
                                 const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out = x;
  out.col(1) = inverse(DBL_EPSILON + (xmat * par + s).array()) +
               (y - inverse(DBL_EPSILON + (xmat * est + s).array())) *
                   ((DBL_EPSILON + (xmat * est + s).array()).pow(-2)) *
                   ((DBL_EPSILON + (xmat * par + s).array()).pow(2));
  return out;
}

Eigen::MatrixXd x0_bin_logit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par,
                             const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out = x;
  out.col(1) =
      y + (logit_linkinv(xmat * par + s) - logit_linkinv(xmat * est + s));
  return out;
}

Eigen::MatrixXd x0_bin_probit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                              const Eigen::Ref<const Eigen::VectorXd> &par,
                              const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c_par =
      probit_linkinv(xmat * par + s) * (1.0 - probit_linkinv(xmat * par + s));
  const Eigen::ArrayXd c_est =
      probit_linkinv(xmat * est + s) * (1.0 - probit_linkinv(xmat * est + s));
  const Eigen::ArrayXd phi_par = 0.5 * M_SQRT1_2 * M_2_SQRTPI *
                                 exp(-0.5 * (xmat * par + s).array().square());
  const Eigen::ArrayXd phi_est = 0.5 * M_SQRT1_2 * M_2_SQRTPI *
                                 exp(-0.5 * (xmat * est + s).array().square());
  Eigen::MatrixXd out = x;
  out.col(1) = probit_linkinv(xmat * par + s) +
               (y - probit_linkinv(xmat * est + s)) * c_par *
                   inverse(DBL_EPSILON + c_est) * phi_est * inverse(phi_par);
  return out;
}

Eigen::MatrixXd x0_bin_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par,
                           const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out = x;
  out.col(1) = log_linkinv(xmat * par + s) +
               (y - log_linkinv(xmat * est + s)) *
                   (inverse(DBL_EPSILON + 1.0 - log_linkinv(xmat * est + s))) *
                   (1.0 - log_linkinv(xmat * par + s));
  return out;
}

Eigen::MatrixXd x0_poi_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par,
                           const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out = x;
  out.col(1) = y + (log_linkinv(xmat * par + s) - log_linkinv(xmat * est + s));
  return out;
}

Eigen::MatrixXd x0_poi_identity(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                const Eigen::Ref<const Eigen::VectorXd> &par,
                                const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out = x;
  out.col(1) = (xmat * par + s).array() * y *
               inverse(DBL_EPSILON + (xmat * est + s).array());
  return out;
}

Eigen::MatrixXd x0_poi_sqrt(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par,
                            const Eigen::Ref<const Eigen::VectorXd> &est)
{
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() *
         (2.0 * y * inverse(DBL_EPSILON + (xmat * par + s).array()) -
          2.0 * (xmat * par + s).array());

  Eigen::MatrixXd out = x;
  out.col(1) = (xmat * par + s).array().square() +
               (xmat * par + s).array() *
                   (inverse(DBL_EPSILON + (xmat * est + s).array()) * y -
                    (xmat * est + s).array());
  return out;
}
