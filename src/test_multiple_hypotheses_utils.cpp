#include "test_multiple_hypotheses_utils.h"
#include "helpers.h"
#include <RcppEigen.h>
#include <functional>
#include <map>
#include <string>
#include <vector>

Eigen::RowVectorXd rmvn(const Eigen::Ref<const Eigen::MatrixXd> &sqrt) {
  Eigen::RowVectorXd u(sqrt.cols());
  for (int i = 0; i < sqrt.cols(); ++i) {
    u(i) = R::rnorm(0, 1.0);
  }
  return u * sqrt;
}

Eigen::MatrixXd w_mean(const Eigen::Ref<const Eigen::MatrixXd> &x,
                       const Eigen::Ref<const Eigen::VectorXd> &par) {
  return Eigen::MatrixXd::Identity(x.cols(), x.cols());
}

Eigen::MatrixXd w_lm(const Eigen::Ref<const Eigen::MatrixXd> &x,
                     const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return static_cast<double>(x.rows()) * (xmat.transpose() * xmat).inverse();
}

Eigen::MatrixXd w_gauss_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c = y * log_linkinv(xmat * par + s) -
                           2.0 * log_linkinv(2.0 * (xmat * par + s));
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_gauss_inverse(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c =
      (DBL_EPSILON + (xmat * par + s).array()).pow(-3) *
      (2.0 * y - 3.0 * inverse(DBL_EPSILON + (xmat * par + s).array()));
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_bin_logit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c =
      logit_linkinv(xmat * par + s) * (1.0 - logit_linkinv(xmat * par + s));
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_bin_probit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd o = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd phi = 0.5 * M_SQRT1_2 * M_2_SQRTPI *
                             exp(-0.5 * (xmat * par + o).array().square());
  const Eigen::ArrayXd dphi = -0.5 * M_SQRT1_2 * M_2_SQRTPI *
                              exp(-0.5 * (xmat * par + o).array().square()) *
                              (xmat * par + o).array();
  const Eigen::ArrayXd s = (y - probit_linkinv(xmat * par + o)) * phi;
  const Eigen::ArrayXd t =
      DBL_EPSILON +
      probit_linkinv(xmat * par + o) * (1.0 - probit_linkinv(xmat * par + o));
  const Eigen::ArrayXd ds =
      -phi.square() + (y - probit_linkinv(xmat * par + o)) * dphi;
  const Eigen::ArrayXd dt = phi * (1.0 - 2.0 * probit_linkinv(xmat * par + o));
  Eigen::ArrayXd c = t.pow(-2) * (ds * t - s * dt);
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_bin_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                          const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c =
      (DBL_EPSILON + 1.0 - log_linkinv(xmat * par + s)).pow(-2) *
      log_linkinv(xmat * par + s) * (y - 1.0);
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_poi_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                          const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c = log_linkinv(xmat * par + s);
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_poi_identity(const Eigen::Ref<const Eigen::MatrixXd> &x,
                               const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c = y * (DBL_EPSILON + (xmat * par + s).array()).pow(-2);
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_poi_sqrt(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c =
      2.0 * y * (DBL_EPSILON + (xmat * par + s).array()).pow(-2) + 2.0;
  return static_cast<double>(x.rows()) *
         (xmat.transpose() * (xmat.array().colwise() * c).matrix()).inverse();
}

Eigen::MatrixXd w_qpoi_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par) {
  const int p = x.cols() - 2;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c =
      -std::pow(phi, -2) * (2.0 * (y - log_linkinv(xmat * beta + s)) +
                            (log_linkinv(-xmat * beta - s)) *
                                square(y - log_linkinv(xmat * beta + s)));
  Eigen::MatrixXd out(p + 1, p + 1);
  out.topLeftCorner(p, p) =
      -(xmat.transpose() * (xmat.array().colwise() *
                            (std::pow(phi, -1) * log_linkinv(xmat * beta + s)))
                               .matrix());
  // out.topRightCorner(p, 1) = Eigen::VectorXd::Zero(p);
  out.topRightCorner(p, 1) =
      -std::pow(phi, -2) *
      (xmat.array().colwise() * (y - log_linkinv(xmat * beta + s)))
          .colwise()
          .sum()
          .transpose();
  out.bottomLeftCorner(1, p) = (xmat.array().colwise() * c).colwise().sum();
  out(p, p) = ((-2.0 * std::pow(phi, -3) * log_linkinv(-xmat * beta - s) *
                square(y - log_linkinv(xmat * beta + s))) +
               std::pow(phi, -2))
                  .sum();
  return static_cast<double>(x.rows()) * out.inverse();
}

Eigen::MatrixXd w_qpoi_identity(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                const Eigen::Ref<const Eigen::VectorXd> &par) {
  const int p = x.cols() - 2;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out(p + 1, p + 1);
  out.topLeftCorner(p, p) = -(
      xmat.transpose() * (xmat.array().colwise() *
                          (std::pow(phi, -1) * y *
                           ((DBL_EPSILON + (xmat * par + s).array()).pow(-2))))
                             .matrix());
  out.topRightCorner(p, 1) =
      -std::pow(phi, -2) *
      (xmat.array().colwise() *
       (y * inverse(DBL_EPSILON + (xmat * par + s).array()) - 1.0))
          .colwise()
          .sum()
          .transpose();
  out.bottomLeftCorner(1, p) =
      -std::pow(phi, -2) *
      (xmat.array().colwise() *
       (square(y) * (DBL_EPSILON + (xmat * par + s).array()).pow(-2) - 1.0))
          .colwise()
          .sum();
  out(p, p) = ((-2.0 * std::pow(phi, -3) *
                inverse(DBL_EPSILON + (xmat * par + s).array()) *
                square(y - (xmat * beta + s).array())) +
               std::pow(phi, -2))
                  .sum();
  return static_cast<double>(x.rows()) * out.inverse();
}

Eigen::MatrixXd w_qpoi_sqrt(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par) {
  const int p = x.cols() - 2;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out(p + 1, p + 1);
  out.topLeftCorner(p, p) =
      -2.0 * std::pow(phi, -1) *
      (xmat.transpose() *
       (xmat.array().colwise() *
        (y * (DBL_EPSILON + (xmat * par + s).array()).pow(-2) + 1.0))
           .matrix());
  out.topRightCorner(p, 1) = -2.0 * std::pow(phi, -2) *
                             (xmat.array().colwise() *
                              (inverse(DBL_EPSILON + (xmat * par + s).array()) *
                               (y - square(y - (xmat * beta + s).array()))))
                                 .colwise()
                                 .sum()
                                 .transpose();
  out.bottomLeftCorner(1, p) =
      -2.0 * std::pow(phi, -2) *
      (xmat.array().colwise() *
       (inverse(DBL_EPSILON + (xmat * par + s).array()) *
        (y - square(y - (xmat * beta + s).array())) *
        ((DBL_EPSILON + (xmat * par + s).array()).pow(-2) *
             (y - square(y - (xmat * beta + s).array())) +
         2.0)))
          .colwise()
          .sum();
  out(p, p) = (-2.0 * std::pow(phi, -3) *
                   (DBL_EPSILON + (xmat * par + s).array()).pow(-2) *
                   square(y - square(y - (xmat * beta + s).array())) +
               std::pow(phi, -2))
                  .sum();
  return static_cast<double>(x.rows()) * out.inverse();
}

Eigen::MatrixXd dg0_inv(const std::string method,
                        const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::VectorXd> &par) {
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
             {"poisson_sqrt", w_poi_sqrt},
             {"quasipoisson_log", w_qpoi_log},
             {"quasipoisson_identity", w_qpoi_identity},
             {"quasipoisson_sqrt", w_qpoi_sqrt}}};
  return w_map[method](x, par);
}

Eigen::MatrixXd shat(const std::string method,
                     const Eigen::Ref<const Eigen::VectorXd> &par,
                     const Eigen::Ref<const Eigen::MatrixXd> &x) {
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd> &,
      const Eigen::Ref<const Eigen::VectorXd> &)>
      g_fn = set_g_fn(method);
  return (1.0 / x.rows()) * ((g_fn(x, par)).transpose() * g_fn(x, par));
}

Eigen::MatrixXd ahat(const Eigen::Ref<const Eigen::MatrixXd> &j,
                     const Eigen::Ref<const Eigen::MatrixXd> &w,
                     const Eigen::Ref<const Eigen::MatrixXd> &s) {
  return (j * w).transpose() * (((j * w) * s * (j * w).transpose()).inverse()) *
         (j * w);
}
