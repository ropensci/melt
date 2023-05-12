#include "helpers.h"
#include <RcppEigen.h>
#include <cfloat>
#include <cmath>
#include <functional>
#include <map>
#include <string>

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &)>
set_g_fn(const std::string method) {
  std::map<std::string, std::function<Eigen::MatrixXd(
                            const Eigen::Ref<const Eigen::MatrixXd> &,
                            const Eigen::Ref<const Eigen::VectorXd> &)>>
      g_map{{{"mean", g_mean},
             {"sd", g_sd},
             {"lm", g_lm},
             {"gaussian_identity", g_lm},
             {"gaussian_log", g_gauss_log},
             {"gaussian_inverse", g_gauss_inverse},
             {"binomial_logit", g_bin_logit},
             {"binomial_probit", g_bin_probit},
             {"binomial_log", g_bin_log},
             {"poisson_log", g_poi_log},
             {"poisson_identity", g_poi_identity},
             {"poisson_sqrt", g_poi_sqrt},
             {"quasipoisson_log", g_qpoi_log},
             {"quasipoisson_identity", g_qpoi_identity},
             {"quasipoisson_sqrt", g_qpoi_sqrt}}};
  return g_map[method];
}

Eigen::VectorXd proj(const Eigen::Ref<const Eigen::MatrixXd> &l,
                     const Eigen::Ref<const Eigen::VectorXd> &x) {
  return x - l.transpose() * ((l * l.transpose()).householderQr().solve(l * x));
}

Eigen::ArrayXd log_linkinv(const Eigen::Ref<const Eigen::VectorXd> &x) {
  return exp(x.array());
}

Eigen::ArrayXd logit_linkinv(const Eigen::Ref<const Eigen::VectorXd> &x) {
  return inverse(1.0 + exp(-x.array()));
}

Eigen::ArrayXd probit_linkinv(const Eigen::Ref<const Eigen::VectorXd> &x) {
  Eigen::ArrayXd out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    out[i] = 0.5 * std::erfc(-x[i] * M_SQRT1_2);
  }
  return out;
}

Eigen::VectorXd mele_mean(const Eigen::Ref<const Eigen::MatrixXd> &x,
                          const Eigen::Ref<const Eigen::ArrayXd> &w) {
  if (w.size() == 0) {
    return x.colwise().mean();
  } else {
    return (w.matrix().transpose() * x) / x.rows();
  }
}

Eigen::VectorXd mele_sd(const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::ArrayXd> &w) {
  if (w.size() == 0) {
    return (x.colwise().mean()).array().sqrt();
  } else {
    return ((w.matrix().transpose() * x) / x.rows()).array().sqrt();
  }
}

Eigen::VectorXd mele_lm(const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::ArrayXd> &w) {
  const Eigen::VectorXd y = x.col(1) - x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  if (w.size() == 0) {
    return xmat.colPivHouseholderQr().solve(y);
  } else {
    const Eigen::MatrixXd wsqrt =
        Eigen::MatrixXd(w.sqrt().matrix().asDiagonal());
    return (wsqrt * xmat).colPivHouseholderQr().solve(wsqrt * y);
  }
}

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::MatrixXd> &x,
                       const Eigen::Ref<const Eigen::VectorXd> &par) {
  return x.rowwise() - par.transpose();
}

Eigen::VectorXd gr_nloglr_mean(const Eigen::Ref<const Eigen::VectorXd> &l,
                               const Eigen::Ref<const Eigen::MatrixXd> &g,
                               const Eigen::Ref<const Eigen::MatrixXd> &x,
                               const Eigen::Ref<const Eigen::VectorXd> &par,
                               const Eigen::Ref<const Eigen::ArrayXd> &w,
                               const bool weighted) {
  Eigen::ArrayXd c(x.rows());
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  }
  return -c.sum() * l;
}

Eigen::MatrixXd g_sd(const Eigen::Ref<const Eigen::MatrixXd> &x,
                     const Eigen::Ref<const Eigen::VectorXd> &par) {
  return x.rowwise() - par.array().square().matrix().transpose();
}

Eigen::VectorXd gr_nloglr_sd(const Eigen::Ref<const Eigen::VectorXd> &l,
                             const Eigen::Ref<const Eigen::MatrixXd> &g,
                             const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par,
                             const Eigen::Ref<const Eigen::ArrayXd> &w,
                             const bool weighted) {
  Eigen::ArrayXd c(x.rows());
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  }
  return -2 * c.sum() * (par.array() * l.array());
}

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::MatrixXd> &x,
                     const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() * (y - (xmat * par + s).array());
}

Eigen::VectorXd gr_nloglr_lm(const Eigen::Ref<const Eigen::VectorXd> &l,
                             const Eigen::Ref<const Eigen::MatrixXd> &g,
                             const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par,
                             const Eigen::Ref<const Eigen::ArrayXd> &w,
                             const bool weighted) {
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c(x.rows());
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  }
  return -(xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_gauss_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() *
         ((y - log_linkinv(xmat * par + s)) * log_linkinv(xmat * par + s));
}

Eigen::VectorXd
gr_nloglr_gauss_log(const Eigen::Ref<const Eigen::VectorXd> &l,
                    const Eigen::Ref<const Eigen::MatrixXd> &g,
                    const Eigen::Ref<const Eigen::MatrixXd> &x,
                    const Eigen::Ref<const Eigen::VectorXd> &par,
                    const Eigen::Ref<const Eigen::ArrayXd> &w,
                    const bool weighted) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c = y * log_linkinv(xmat * par + s) -
                     2.0 * log_linkinv(2.0 * (xmat * par + s));
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_gauss_inverse(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() *
         (-(y - inverse(DBL_EPSILON + (xmat * par + s).array())) *
          (DBL_EPSILON + (xmat * par + s).array()).pow(-2));
}

Eigen::VectorXd
gr_nloglr_gauss_inverse(const Eigen::Ref<const Eigen::VectorXd> &l,
                        const Eigen::Ref<const Eigen::MatrixXd> &g,
                        const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::VectorXd> &par,
                        const Eigen::Ref<const Eigen::ArrayXd> &w,
                        const bool weighted) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c =
      (DBL_EPSILON + (xmat * par + s).array()).pow(-3) *
      (2.0 * y - 3.0 * inverse(DBL_EPSILON + (xmat * par + s).array()));
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_bin_logit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() * (y - logit_linkinv(xmat * par + s));
}

Eigen::VectorXd
gr_nloglr_bin_logit(const Eigen::Ref<const Eigen::VectorXd> &l,
                    const Eigen::Ref<const Eigen::MatrixXd> &g,
                    const Eigen::Ref<const Eigen::MatrixXd> &x,
                    const Eigen::Ref<const Eigen::VectorXd> &par,
                    const Eigen::Ref<const Eigen::ArrayXd> &w,
                    const bool weighted) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c =
      -logit_linkinv(xmat * par + s) * (1.0 - logit_linkinv(xmat * par + s));
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_bin_probit(const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  const Eigen::ArrayXd c =
      inverse(DBL_EPSILON + probit_linkinv(xmat * par + s) *
                                (1.0 - probit_linkinv(xmat * par + s)));
  const Eigen::ArrayXd phi = 0.5 * M_SQRT1_2 * M_2_SQRTPI *
                             exp(-0.5 * (xmat * par + s).array().square());
  return xmat.array().colwise() *
         (c * (y - probit_linkinv(xmat * par + s)) * phi);
}

Eigen::VectorXd
gr_nloglr_bin_probit(const Eigen::Ref<const Eigen::VectorXd> &l,
                     const Eigen::Ref<const Eigen::MatrixXd> &g,
                     const Eigen::Ref<const Eigen::MatrixXd> &x,
                     const Eigen::Ref<const Eigen::VectorXd> &par,
                     const Eigen::Ref<const Eigen::ArrayXd> &w,
                     const bool weighted) {
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
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_bin_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                          const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() *
         ((inverse(DBL_EPSILON + 1.0 - log_linkinv(xmat * par + s))) *
          (y - log_linkinv(xmat * par + s)));
}

Eigen::VectorXd gr_nloglr_bin_log(const Eigen::Ref<const Eigen::VectorXd> &l,
                                  const Eigen::Ref<const Eigen::MatrixXd> &g,
                                  const Eigen::Ref<const Eigen::MatrixXd> &x,
                                  const Eigen::Ref<const Eigen::VectorXd> &par,
                                  const Eigen::Ref<const Eigen::ArrayXd> &w,
                                  const bool weighted) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c = (DBL_EPSILON + 1.0 - log_linkinv(xmat * par + s)).pow(-2) *
                     log_linkinv(xmat * par + s) * (y - 1.0);
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_poi_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                          const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() * (y - log_linkinv(xmat * par + s));
}

Eigen::VectorXd gr_nloglr_poi_log(const Eigen::Ref<const Eigen::VectorXd> &l,
                                  const Eigen::Ref<const Eigen::MatrixXd> &g,
                                  const Eigen::Ref<const Eigen::MatrixXd> &x,
                                  const Eigen::Ref<const Eigen::VectorXd> &par,
                                  const Eigen::Ref<const Eigen::ArrayXd> &w,
                                  const bool weighted) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c = -log_linkinv(xmat * par + s);
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_poi_identity(const Eigen::Ref<const Eigen::MatrixXd> &x,
                               const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() *
         (y * inverse(DBL_EPSILON + (xmat * par + s).array()) - 1.0);
}

Eigen::VectorXd
gr_nloglr_poi_identity(const Eigen::Ref<const Eigen::VectorXd> &l,
                       const Eigen::Ref<const Eigen::MatrixXd> &g,
                       const Eigen::Ref<const Eigen::MatrixXd> &x,
                       const Eigen::Ref<const Eigen::VectorXd> &par,
                       const Eigen::Ref<const Eigen::ArrayXd> &w,
                       const bool weighted) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c = -y * (DBL_EPSILON + (xmat * par + s).array()).pow(-2);
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_poi_sqrt(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  return xmat.array().colwise() *
         (2.0 * y * inverse(DBL_EPSILON + (xmat * par + s).array()) -
          2.0 * (xmat * par + s).array());
}

Eigen::VectorXd gr_nloglr_poi_sqrt(const Eigen::Ref<const Eigen::VectorXd> &l,
                                   const Eigen::Ref<const Eigen::MatrixXd> &g,
                                   const Eigen::Ref<const Eigen::MatrixXd> &x,
                                   const Eigen::Ref<const Eigen::VectorXd> &par,
                                   const Eigen::Ref<const Eigen::ArrayXd> &w,
                                   const bool weighted) {
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c =
      -2.0 * y * (DBL_EPSILON + (xmat * par + s).array()).pow(-2) - 2.0;
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_qpoi_log(const Eigen::Ref<const Eigen::MatrixXd> &x,
                           const Eigen::Ref<const Eigen::VectorXd> &par) {
  const int p = x.cols() - 2;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out(x.rows(), p + 1);
  out.leftCols(p) = (1.0 / phi) * (xmat.array().colwise() *
                                   (y - log_linkinv(xmat * beta + s)));
  out.col(p) = inverse(phi * phi * log_linkinv(xmat * beta + s)) *
                   square(y - log_linkinv(xmat * beta + s)) -
               1.0 / phi;
  return out;
}

Eigen::VectorXd gr_nloglr_qpoi_log(const Eigen::Ref<const Eigen::VectorXd> &l,
                                   const Eigen::Ref<const Eigen::MatrixXd> &g,
                                   const Eigen::Ref<const Eigen::MatrixXd> &x,
                                   const Eigen::Ref<const Eigen::VectorXd> &par,
                                   const Eigen::Ref<const Eigen::ArrayXd> &w,
                                   const bool weighted) {
  const int p = x.cols() - 2;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c(x.rows());
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  }
  const Eigen::ArrayXd c2 =
      -std::pow(phi, -2) * (c * (2.0 * (y - log_linkinv(xmat * beta + s)) +
                                 (log_linkinv(-xmat * beta - s)) *
                                     square(y - log_linkinv(xmat * beta + s))));
  Eigen::MatrixXd out(p + 1, p + 1);
  out.topLeftCorner(p, p) =
      -(xmat.transpose() *
        (xmat.array().colwise() *
         (std::pow(phi, -1) * c * log_linkinv(xmat * beta + s)))
            .matrix());
  out.topRightCorner(p, 1) =
      -std::pow(phi, -2) *
      (xmat.array().colwise() * (c * (y - log_linkinv(xmat * beta + s))))
          .colwise()
          .sum()
          .transpose();
  out.bottomLeftCorner(1, p) = (xmat.array().colwise() * c2).colwise().sum();
  out(p, p) = (c * ((-2.0 * std::pow(phi, -3) * log_linkinv(-xmat * beta - s) *
                     square(y - log_linkinv(xmat * beta + s))) +
                    std::pow(phi, -2)))
                  .sum();
  return out * l;
}

Eigen::MatrixXd g_qpoi_identity(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                const Eigen::Ref<const Eigen::VectorXd> &par) {
  const int p = x.cols() - 2;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out(x.rows(), p + 1);
  out.leftCols(p) =
      (1.0 / phi) *
      (xmat.array().colwise() *
       (y * inverse(DBL_EPSILON + (xmat * par + s).array()) - 1.0));
  out.col(p) = std::pow(phi, -2) *
                   inverse(DBL_EPSILON + (xmat * par + s).array()) *
                   square(y - (xmat * beta + s).array()) -
               1.0 / phi;
  return out;
}

Eigen::VectorXd
gr_nloglr_qpoi_identity(const Eigen::Ref<const Eigen::VectorXd> &l,
                        const Eigen::Ref<const Eigen::MatrixXd> &g,
                        const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::VectorXd> &par,
                        const Eigen::Ref<const Eigen::ArrayXd> &w,
                        const bool weighted) {
  const int p = x.cols() - 2;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c(x.rows());
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  }
  Eigen::MatrixXd out(p + 1, p + 1);
  out.topLeftCorner(p, p) = -(
      xmat.transpose() * (xmat.array().colwise() *
                          (std::pow(phi, -1) * c * y *
                           ((DBL_EPSILON + (xmat * par + s).array()).pow(-2))))
                             .matrix());
  out.topRightCorner(p, 1) =
      -std::pow(phi, -2) *
      (xmat.array().colwise() *
       (c * (y * inverse(DBL_EPSILON + (xmat * par + s).array()) - 1.0)))
          .colwise()
          .sum()
          .transpose();
  out.bottomLeftCorner(1, p) =
      -std::pow(phi, -2) *
      (xmat.array().colwise() *
       (c *
        (square(y) * (DBL_EPSILON + (xmat * par + s).array()).pow(-2) - 1.0)))
          .colwise()
          .sum();
  out(p, p) = (c * ((-2.0 * std::pow(phi, -3) *
                     inverse(DBL_EPSILON + (xmat * par + s).array()) *
                     square(y - (xmat * beta + s).array())) +
                    std::pow(phi, -2)))
                  .sum();
  return out * l;
}

Eigen::MatrixXd g_qpoi_sqrt(const Eigen::Ref<const Eigen::MatrixXd> &x,
                            const Eigen::Ref<const Eigen::VectorXd> &par) {
  // The first two columns of x are the offset and y
  const int p = x.cols() - 2;
  // The last element of par is the dispersion
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::MatrixXd out(x.rows(), p + 1);
  out.leftCols(p) =
      xmat.array().colwise() *
      ((2.0 / phi) * (y * inverse(DBL_EPSILON + (xmat * par + s).array()) -
                      (xmat * par + s).array()));
  out.col(p) = std::pow(phi, -2) *
                   ((DBL_EPSILON + (xmat * par + s).array()).pow(-2)) *
                   square(y - (xmat * beta + s).array().pow(2)) -
               1.0 / phi;
  return out;
}

Eigen::VectorXd
gr_nloglr_qpoi_sqrt(const Eigen::Ref<const Eigen::VectorXd> &l,
                    const Eigen::Ref<const Eigen::MatrixXd> &g,
                    const Eigen::Ref<const Eigen::MatrixXd> &x,
                    const Eigen::Ref<const Eigen::VectorXd> &par,
                    const Eigen::Ref<const Eigen::ArrayXd> &w,
                    const bool weighted) {
  // The first two columns of x are the offset and y
  const int p = x.cols() - 2;
  // The last element of par is the dispersion
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::VectorXd s = x.col(0);
  const Eigen::ArrayXd y = x.col(1);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 2);
  Eigen::ArrayXd c(x.rows());
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  }
  Eigen::MatrixXd out(p + 1, p + 1);
  out.topLeftCorner(p, p) =
      -2.0 * std::pow(phi, -1) *
      (xmat.transpose() *
       (xmat.array().colwise() *
        (c * (y * (DBL_EPSILON + (xmat * par + s).array()).pow(-2) + 1.0)))
           .matrix());
  out.topRightCorner(p, 1) =
      -2.0 * std::pow(phi, -2) *
      (xmat.array().colwise() *
       (c * inverse(DBL_EPSILON + (xmat * par + s).array()) *
        (y - square(y - (xmat * beta + s).array()))))
          .colwise()
          .sum()
          .transpose();
  out.bottomLeftCorner(1, p) =
      -2.0 * std::pow(phi, -2) *
      (xmat.array().colwise() *
       (c * inverse(DBL_EPSILON + (xmat * par + s).array()) *
        (y - square(y - (xmat * beta + s).array())) *
        ((DBL_EPSILON + (xmat * par + s).array()).pow(-2) *
             (y - square(y - (xmat * beta + s).array())) +
         2.0)))
          .colwise()
          .sum();
  out(p, p) = (c * ((-2.0 * std::pow(phi, -3) *
                     (DBL_EPSILON + (xmat * par + s).array()).pow(-2) *
                     square(y - square(y - (xmat * beta + s).array()))) +
                    std::pow(phi, -2)))
                  .sum();
  return out * l;
}
