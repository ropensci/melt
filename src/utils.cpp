#include "utils.h"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <functional>
#include <map>
#include <string>

// [[Rcpp::export]]
int getMaxThreads() {
  #ifdef _OPENMP
  return omp_get_max_threads();
  #else
  return 1;
  #endif
}

// [[Rcpp::export]]
int getRank(const Eigen::Map<Eigen::MatrixXd> &x) {
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  return lu_decomp.rank();
}

Eigen::ArrayXd inverse_linkinv(const Eigen::Ref<const Eigen::VectorXd> &x) {
  return x.array().inverse();
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
};
Eigen::VectorXd mele_lm(const Eigen::Ref<const Eigen::MatrixXd> &data,
                        const Eigen::Ref<const Eigen::ArrayXd> &w) {
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
  return c.sum() * l;
}

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::MatrixXd> &x,
                     const Eigen::Ref<const Eigen::VectorXd> &par) {
  return x.rightCols(x.cols() - 1).array().colwise() *
         (x.col(0) - x.rightCols(x.cols() - 1) * par).array();
}
Eigen::VectorXd gr_nloglr_lm(const Eigen::Ref<const Eigen::VectorXd> &l,
                             const Eigen::Ref<const Eigen::MatrixXd> &g,
                             const Eigen::Ref<const Eigen::MatrixXd> &x,
                             const Eigen::Ref<const Eigen::VectorXd> &par,
                             const Eigen::Ref<const Eigen::ArrayXd> &w,
                             const bool weighted){
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
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
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  return xmat.array().colwise() *
         ((y - log_linkinv(xmat * par)) * log_linkinv(xmat * par));
}
Eigen::VectorXd gr_nloglr_gauss_log(
    const Eigen::Ref<const Eigen::VectorXd> &l,
    const Eigen::Ref<const Eigen::MatrixXd> &g,
    const Eigen::Ref<const Eigen::MatrixXd> &x,
    const Eigen::Ref<const Eigen::VectorXd> &par,
    const Eigen::Ref<const Eigen::ArrayXd> &w,
    const bool weighted) {
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  Eigen::ArrayXd c = (y * log_linkinv(xmat * par) -
                      2.0 * log_linkinv(2.0 * xmat * par)) *
                     inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  if (weighted){
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_gauss_inverse(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                const Eigen::Ref<const Eigen::VectorXd> &par) {
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  return xmat.array().colwise() *
         (-(y - inverse_linkinv(xmat * par)) *
         inverse_linkinv(xmat * par).square());
}
Eigen::VectorXd gr_nloglr_gauss_inverse(
    const Eigen::Ref<const Eigen::VectorXd> &l,
    const Eigen::Ref<const Eigen::MatrixXd> &g,
    const Eigen::Ref<const Eigen::MatrixXd> &x,
    const Eigen::Ref<const Eigen::VectorXd> &par,
    const Eigen::Ref<const Eigen::ArrayXd> &w,
    const bool weighted) {
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  Eigen::ArrayXd c = inverse_linkinv(xmat * par).cube() *
                     (2.0 * y - 3.0 * inverse_linkinv(xmat * par));
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}





Eigen::MatrixXd g_bin_logit(const Eigen::Ref<const Eigen::MatrixXd>& x,
                            const Eigen::Ref<const Eigen::VectorXd>& par) {
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  return xmat.array().colwise() * (y - logit_linkinv(xmat * par));
}
Eigen::VectorXd gr_nloglr_bin_logit(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w,
    const bool weighted) {
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  Eigen::ArrayXd c = -logit_linkinv(xmat * par) *
    (1.0 - logit_linkinv(xmat * par));
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_bin_probit(const Eigen::Ref<const Eigen::MatrixXd>& x,
                             const Eigen::Ref<const Eigen::VectorXd>& par) {
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd phi = 0.5 * M_SQRT1_2 * M_2_SQRTPI *
    exp(-0.5 * (xmat * par).array().square());
  return xmat.array().colwise() *
    ((y * inverse(probit_linkinv(xmat * par)) - 1.0) * phi);
}
Eigen::VectorXd gr_nloglr_bin_probit(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w,
    const bool weighted) {
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd phi = 0.5 * M_SQRT1_2 * M_2_SQRTPI *
    exp(-0.5 * (xmat * par).array().square());
  const Eigen::ArrayXd dphi = -0.5 * M_SQRT1_2 * M_2_SQRTPI *
    exp(-0.5 * (xmat * par).array().square()) * (xmat * par).array();
  Eigen::ArrayXd c = y * (dphi * inverse(probit_linkinv(xmat * par)) -
    square(inverse(probit_linkinv(xmat * par)) * phi)) - dphi;
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}

Eigen::MatrixXd g_bin_log(const Eigen::Ref<const Eigen::MatrixXd>& x,
                          const Eigen::Ref<const Eigen::VectorXd>& par) {
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  return xmat.array().colwise() * ((inverse(1.0 - log_linkinv(xmat * par))) *
                    (y - log_linkinv(xmat * par)));
}
Eigen::VectorXd gr_nloglr_bin_log(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w,
    const bool weighted) {
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  Eigen::ArrayXd c = square((1.0 - log_linkinv(xmat * par)).inverse()) *
    log_linkinv(xmat * par) * (y - 1.0);
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * c;
  }
  return (xmat.transpose() * (xmat.array().colwise() * c).matrix()) * l;
}


















std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &)>
g_fn2(const std::string method) {
  std::map<std::string, std::function<Eigen::MatrixXd(
                            const Eigen::Ref<const Eigen::MatrixXd> &,
                            const Eigen::Ref<const Eigen::VectorXd> &)>>
      g_map{{{"mean", g_mean},
             {"lm", g_lm},
             {"gaussian_identity", g_lm},
             {"gaussian_log", g_gauss_log},
             {"gaussian_inverse", g_gauss_inverse},
             {"binomial_logit", g_bin_logit},
             {"binomial_probit", g_bin_probit},
             {"binomial_log", g_bin_log}}};
  return g_map[method];
}

double computeQuantile(const Rcpp::NumericVector &x, const double prob) {
  Rcpp::Environment stats("package:stats");
  Rcpp::Function quantile = stats["quantile"];
  return Rcpp::as<double>(quantile(x, Rcpp::Named("probs") = prob));
}
