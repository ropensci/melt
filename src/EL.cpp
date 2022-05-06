#include "EL.h"

// lm
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
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w,
    const bool weighted)
{
  const double n = static_cast<double>(g.rows());
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd denom = Eigen::VectorXd::Ones(g.rows()) + g * l;
  if (weighted) {
    const Eigen::MatrixXd xx = x.array().colwise() * (w / denom);
    return -(x.transpose() * xx) * l;
  } else {
    const Eigen::MatrixXd xx = x.array().colwise() / denom;
    return -(x.transpose() * xx) * l;
  }
}


// binomial family
Eigen::MatrixXd g_bin_logit(const Eigen::Ref<const Eigen::MatrixXd>& x,
                            const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  return xmat.array().colwise() * (y - logit_linkinv(xmat * par));
}
Eigen::VectorXd gr_nloglr_bin_logit(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w,
    const bool weighted)
{
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd num =
    logit_linkinv(x * par) * (1.0 - logit_linkinv(x * par));
  const Eigen::ArrayXd denom = Eigen::VectorXd::Ones(g.rows()) + g * l;
  if (weighted) {
    const Eigen::MatrixXd cx = x.array().colwise() * (w * num / denom);
    return -(x.transpose() * cx) * l;
  } else {
    const Eigen::MatrixXd cx = x.array().colwise() * (num / denom);
    return -(x.transpose() * cx) * l;
  }
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
    const Eigen::Ref<const Eigen::ArrayXd>& w,
    const bool weighted)
{
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd num =
    -exp(-(x * par).array().square() * 0.5) * M_SQRT1_2 * M_2_SQRTPI * 0.5;
  const Eigen::ArrayXd denom = Eigen::VectorXd::Ones(g.rows()) + g * l;
  if (weighted) {
    const Eigen::MatrixXd cx = x.array().colwise() * (w * num / denom);
    return -(x.transpose() * cx) * l;
  } else {
    const Eigen::MatrixXd cx = x.array().colwise() * (num / denom);
    return -(x.transpose() * cx) * l;
  }
}

Eigen::MatrixXd g_bin_log(const Eigen::Ref<const Eigen::MatrixXd>& x,
                          const Eigen::Ref<const Eigen::VectorXd>& par)
{
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
    const bool weighted)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c = square((1.0 - log_linkinv(xmat * par)).inverse()) *
    log_linkinv(xmat * par) * (y - 1.0) *
    inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  if (weighted) {
    const Eigen::MatrixXd cx = xmat.array().colwise() * (w * c);
    return (xmat.transpose() * cx) * l;
  } else {
    const Eigen::MatrixXd cx = xmat.array().colwise() * c;
    return (xmat.transpose() * cx) * l;
  }
}


// poisson family
Eigen::MatrixXd g_poi_log(const Eigen::Ref<const Eigen::MatrixXd>& data,
                          const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const Eigen::ArrayXd y = data.col(0);
  const Eigen::MatrixXd xmat = data.rightCols(data.cols() - 1);
  return xmat.array().colwise() * (y - log_linkinv(xmat * par));
}
Eigen::VectorXd gr_nloglr_poi_log(const Eigen::Ref<const Eigen::VectorXd>& l,
                                  const Eigen::Ref<const Eigen::MatrixXd>& g,
                                  const Eigen::Ref<const Eigen::MatrixXd>& data,
                                  const Eigen::Ref<const Eigen::VectorXd>& par,
                                  const Eigen::Ref<const Eigen::ArrayXd>& w,
                                  const bool weighted)
{
  const Eigen::MatrixXd xmat = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd num = log_linkinv(xmat * par);
  const Eigen::ArrayXd denom = Eigen::VectorXd::Ones(g.rows()) + (g * l);
  if (weighted) {
    const Eigen::MatrixXd cx = xmat.array().colwise() * (w * num / denom);
    return -(xmat.transpose() * cx) * l;
  } else {
    const Eigen::MatrixXd cx = xmat.array().colwise() * (num / denom);
    return -(xmat.transpose() * cx) * l;
  }
}

Eigen::MatrixXd g_poi_identity(const Eigen::Ref<const Eigen::MatrixXd>& x,
                               const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  return xmat.array().colwise() *
    (inverse((xmat * par).array()) * (y - (xmat * par).array()));
}
Eigen::VectorXd gr_nloglr_poi_identity(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w,
    const bool weighted)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c =
    -inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) * y *
    square(((xmat * par).array().inverse()));
  if (weighted) {
    const Eigen::MatrixXd cx = xmat.array().colwise() * (w * c);
    return (xmat.transpose() * cx) * l;
  } else {
    const Eigen::MatrixXd cx = xmat.array().colwise() * c;
    return (xmat.transpose() * cx) * l;
  }
}

Eigen::MatrixXd g_poi_sqrt(const Eigen::Ref<const Eigen::MatrixXd>& x,
                           const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  return xmat.array().colwise() *
    (2.0 * inverse((xmat * par).array()) * (y - square((xmat * par).array())));
}
Eigen::VectorXd gr_nloglr_poi_sqrt(const Eigen::Ref<const Eigen::VectorXd>& l,
                                   const Eigen::Ref<const Eigen::MatrixXd>& g,
                                   const Eigen::Ref<const Eigen::MatrixXd>& x,
                                   const Eigen::Ref<const Eigen::VectorXd>& par,
                                   const Eigen::Ref<const Eigen::ArrayXd>& w,
                                   const bool weighted)
{
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  const Eigen::ArrayXd c =
    -2.0 * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array()) *
    (y * square(((xmat * par).array().inverse())) + 1.0);
  if (weighted) {
    const Eigen::MatrixXd cx = xmat.array().colwise() * (w * c);
    return (xmat.transpose() * cx) * l;
  } else {
    const Eigen::MatrixXd cx = xmat.array().colwise() * c;
    return (xmat.transpose() * cx) * l;
  }
}


// quasibinomial family
Eigen::MatrixXd g_qbin_logit(const Eigen::Ref<const Eigen::MatrixXd>& x,
                             const Eigen::Ref<const Eigen::VectorXd>& par)
{
  const int p = x.cols() - 1;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(p);

  Eigen::MatrixXd out(x.rows(), p + 1);
  out.leftCols(p) = xmat.array().colwise() * (y - logit_linkinv(xmat * beta));
  // out.rightCols(1) = square(y - logit_linkinv(xmat * beta)) *
  //   inverse(phi * phi * logit_linkinv(xmat * beta) *
  //   (1.0 - logit_linkinv(xmat * beta))) - 1.0 / phi;
  out.col(p) = inverse(phi * phi * logit_linkinv(xmat * beta) *
    (1.0 - logit_linkinv(xmat * beta))) *
    square(y - logit_linkinv(xmat * beta)) - 1.0 / phi;
  return out;
}
Eigen::VectorXd gr_nloglr_qbin_logit(
    const Eigen::Ref<const Eigen::VectorXd>& l,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::VectorXd>& par,
    const Eigen::Ref<const Eigen::ArrayXd>& w,
    const bool weighted)
{
  const int p = x.cols() - 1;
  const Eigen::VectorXd beta = par.head(p);
  const double phi = par(p);
  const Eigen::ArrayXd y = x.col(0);
  const Eigen::MatrixXd xmat = x.rightCols(x.cols() - 1);
  Eigen::ArrayXd c(x.rows());
  if (weighted) {
    c = w * inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  } else {
    c = inverse((Eigen::VectorXd::Ones(g.rows()) + g * l).array());
  }

  Eigen::ArrayXd c2 = inverse(phi * phi * logit_linkinv(xmat * beta) *
    (1.0 - logit_linkinv(xmat * beta))) *
    (square(logit_linkinv(xmat * beta)) * (1.0 - 2.0 * y) +
    square(y) * (2.0 * logit_linkinv(xmat * beta) - 1.0));

  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(p + 1, p + 1);
  out.topLeftCorner(p, p) = xmat.transpose() *
    (xmat.array().colwise() * c).matrix();
  out.bottomLeftCorner(1, p) = (xmat.array().colwise() * c2).colwise().sum();
  out.topRightCorner(p, 1) =
    (xmat.array().colwise() * c2).colwise().sum().transpose();
  out(p, p) = (c * (-2.0 * logit_linkinv(xmat * beta) *
    (1.0 - logit_linkinv(xmat * beta)) *
    square(y - logit_linkinv(xmat * beta)) * std::pow(phi, -3) +
    std::pow(phi, -2))).sum();
  return out * l;
}


/* EL class (evaluation)
 * Last updated: 04/07/22
 */
EL::EL(
  const std::string method,
  const Eigen::Ref<const Eigen::VectorXd>& par0,
  const Eigen::Ref<const Eigen::MatrixXd>& x,
  const int maxit_l,
  const double tol_l,
  const double th,
  const Eigen::Ref<const Eigen::ArrayXd>& wt)
  : par{par0},
    l{Eigen::VectorXd::Zero(par0.size())},
    mele_fn{set_mele_fn(method)},
    w{wt},
    maxit_l{maxit_l},
    tol_l{tol_l},
    th{th},
    n{static_cast<int>(x.rows())},
    g_fn{EL::set_g_fn(method)}
{
  set_el(g_fn(x, par), wt);
}
EL::EL(
  const Eigen::Ref<const Eigen::MatrixXd>& g,
  const int maxit_l,
  const double tol_l,
  const double th,
  const Eigen::Ref<const Eigen::ArrayXd>& wt)
  : par{},
    l{Eigen::VectorXd::Zero(g.cols())},
    mele_fn{},
    w{wt},
    maxit_l{maxit_l},
    tol_l{tol_l},
    th{th},
    n{static_cast<int>(g.rows())},
    g_fn{}
{
  set_el(g, wt);
}

std::function<Eigen::VectorXd(const Eigen::Ref<const Eigen::MatrixXd>&,
                              const Eigen::Ref<const Eigen::ArrayXd>&)>
  EL::set_mele_fn(const std::string method)
{
  std::map<std::string, std::function<Eigen::VectorXd(
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::ArrayXd>&)>>
        mele_map{{{"mean", mele_mean},
                  {"lm", mele_lm},
                  {"gaussian_identity", mele_lm},
                  {"gaussian_log", mele_lm},
                  {"gaussian_inverse", mele_lm},
                  {"binomial_logit", mele_lm},
                  {"binomial_probit", mele_lm},
                  {"binomial_log", mele_lm},
                  {"poisson_log", mele_lm},
                  {"poisson_identity", mele_lm},
                  {"poisson_sqrt", mele_lm},
                  {"quasibinomial_logit", mele_lm}}};
  return mele_map[method];
}

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd>&,
                              const Eigen::Ref<const Eigen::VectorXd>&)>
  EL::set_g_fn(const std::string method)
{
  std::map<std::string, std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&)>>
        g_map{{{"mean", g_mean},
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
               {"quasibinomial_logit", g_qbin_logit}}};
  return g_map[method];
}

void EL::set_el(const Eigen::Ref<const Eigen::MatrixXd>& g,
                const Eigen::Ref<const Eigen::ArrayXd>& w)
{
  // maximization
  while (!conv && iter != maxit_l && nllr <= th) {
    // pseudo log
    const PSEUDO_LOG pl(Eigen::VectorXd::Ones(n) + g * l, w);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * pl.sqrt_neg_d2plog;
    // propose new lambda by NR method with least square
    Eigen::VectorXd step = (J.transpose() * J).ldlt().solve(
      J.transpose() * (pl.dplog / pl.sqrt_neg_d2plog).matrix());
    // update function value
    nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step), w);
    // step halving to ensure increase in function value
    if (nllr < pl.plog_sum) {
      step /= 2;
      nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step), w);
    }
    // convergence check
    if (step.norm() < tol_l * l.norm() + tol_l * tol_l) {
      conv = true;
    }
    ++iter;
    // update lambda
    l += step;
  }
}

Eigen::ArrayXd EL::logp_g(const Eigen::Ref<const Eigen::MatrixXd>& g) const
{
  if (w.size() == 0) {
    return
    -log(n) - PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g * l);
  } else {
    return log(w) - log(n) -
      PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g * l, w);
  }
}

Eigen::ArrayXd EL::logp(const Eigen::Ref<const Eigen::MatrixXd>& x) const
{
  if (w.size() == 0) {
    return
    -log(n) - PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g_fn(x, par) * l);
  } else {
    return log(w) - log(n) -
      PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g_fn(x, par) * l, w);
  }
}

double EL::loglik() const
{
  if (w.size() == 0) {
    return -nllr - n * log(n);
  } else {
    return -nllr - (w * (log(n) - log(w))).sum();
  }
}


/* MINEL class (minimization)
 * Last updated: 04/20/22
 */
MINEL::MINEL(const std::string method,
             const Eigen::Ref<const Eigen::VectorXd>& par0,
             const Eigen::Ref<const Eigen::MatrixXd>& x,
             const Eigen::Ref<const Eigen::MatrixXd>& lhs,
             const Eigen::Ref<const Eigen::VectorXd>& rhs,
             const int maxit,
             const int maxit_l,
             const double tol,
             const double tol_l,
             const double step,
             const double th,
             const Eigen::Ref<const Eigen::ArrayXd>& wt)
  : maxit{maxit},
    maxit_l{maxit_l},
    tol{tol},
    tol_l{tol_l},
    gamma{step},
    th{th},
    n{static_cast<int>(x.rows())},
    weighted{wt.size() != 0},
    g_fn{MINEL::set_g_fn(method)},
    gr_fn{MINEL::set_gr_fn(method)}
{
  /// initialization ///
  // orthogonal projection matrix
  const Eigen::MatrixXd proj =
    Eigen::MatrixXd::Identity(lhs.cols(), lhs.cols()) -
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * lhs;
  // parameter (constraint imposed)
  par = proj * par0 + lhs.transpose() * (lhs * lhs.transpose()).inverse() * rhs;
  // estimating function
  Eigen::MatrixXd g = g_fn(x, par);
  // lambda
  l = EL(g, maxit_l, tol_l, th, wt).l;
  // function value (-logLR)
  nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * l, wt);
  // function norm
  const double norm0 = (proj * gr_fn(l, g, x, par, wt, weighted)).norm();

  /// minimization (projected gradient descent) ///
  // double gamma = 1.0;
  while (!conv && iter != maxit && nllr <= th) {
    // update parameter
    Eigen::VectorXd par_tmp =
      par - gamma * proj * gr_fn(l, g, x, par, wt, weighted);
    // update estimating function
    Eigen::MatrixXd g_tmp = g_fn(x, par_tmp);
    // update lambda
    Eigen::VectorXd l_tmp = EL(g_tmp, maxit_l, tol_l, th, wt).l;
    // update function value
    const double f0 = nllr;
    nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g_tmp * l_tmp, wt);
    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < nllr || nllr < 0) {
      // reduce step size
      gamma /= 2;
      if (gamma < 1e-30) {
        break;
      }
      // propose new parameter
      par_tmp = par - gamma * proj * gr_fn(l, g, x, par, wt, weighted);
      // propose new lambda
      g_tmp = g_fn(x, par_tmp);
      l_tmp = EL(g_tmp, maxit_l, tol_l, th, wt).l;
      // propose new function value
      nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g_tmp * l_tmp, wt);
    }

    if (nllr >= f0) {
      nllr = f0;
      ++iter;
      break;
    } else if (std::isnan(nllr)) {
      // return initial values
      par =
        proj * par0 + lhs.transpose() * (lhs * lhs.transpose()).inverse() * rhs;
      Eigen::MatrixXd g = g_fn(x, par);
      l = EL(g, maxit_l, tol_l, th, wt).l;
      nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * l, wt);
    }

    // update
    const double step = (par - par_tmp).norm();
    par = std::move(par_tmp);
    l = std::move(l_tmp);
    g = std::move(g_tmp);
    // convergence check
    // if ((proj * gr_fn(l, g, x, par, wt, weighted)).norm() <
    //   tol * norm0 + tol * tol || step < tol * par.norm() + tol * tol) {
    //   conv = true;
    // }
    if ((proj * gr_fn(l, g, x, par, wt, weighted)).norm() < tol ||
        step < tol * par.norm() + tol * tol) {
      conv = true;
    }
    // if ((proj * gr_fn(l, g, x, par, wt, weighted)).norm() < tol) {
    //   conv = true;
    // }
    ++iter;
  }
}





std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd>&,
                              const Eigen::Ref<const Eigen::VectorXd>&)>
  MINEL::set_g_fn(const std::string method)
{
  std::map<std::string, std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&)>>
        g_map{{{"mean", g_mean},
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
               {"quasibinomial_logit", g_qbin_logit}}};
  return g_map[method];
}

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::VectorXd>&,
                              const Eigen::Ref<const Eigen::MatrixXd>&,
                              const Eigen::Ref<const Eigen::MatrixXd>&,
                              const Eigen::Ref<const Eigen::VectorXd>&,
                              const Eigen::Ref<const Eigen::ArrayXd>&,
                              const bool)>
  MINEL::set_gr_fn(const std::string method)
{
  std::map<std::string, std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::VectorXd>&,
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&,
      const Eigen::Ref<const Eigen::ArrayXd>&,
      const bool)>> gr_map{
        {{"mean", gr_nloglr_mean},
         {"lm", gr_nloglr_lm},
         {"gaussian_identity", gr_nloglr_lm},
         {"gaussian_log", gr_nloglr_gauss_log},
         {"gaussian_inverse", gr_nloglr_gauss_inverse},
         {"binomial_logit", gr_nloglr_bin_logit},
         {"binomial_probit", gr_nloglr_bin_probit},
         {"binomial_log", gr_nloglr_bin_log},
         {"poisson_log", gr_nloglr_poi_log},
         {"poisson_identity", gr_nloglr_poi_identity},
         {"poisson_sqrt", gr_nloglr_poi_sqrt},
         {"quasibinomial_logit", gr_nloglr_qbin_logit}}};
  return gr_map[method];
}

Eigen::ArrayXd MINEL::logp(const Eigen::Ref<const Eigen::MatrixXd>& x,
                           const Eigen::Ref<const Eigen::ArrayXd>& wt) const
{
  if (weighted) {
    return log(wt) - log(n) -
      PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g_fn(x, par) * l, wt);
  } else {
    return
    -log(n) - PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g_fn(x, par) * l);
  }
}

double MINEL::loglik(const Eigen::Ref<const Eigen::ArrayXd>& wt) const
{
  if (weighted) {
    return -nllr - (wt * (log(n) - log(wt))).sum();
  } else {
    return -nllr - n * log(n);
  }
}


/* PSEUDO_LOG class
 * Last updated: 04/07/22
 *
 */
PSEUDO_LOG::PSEUDO_LOG(const Eigen::Ref<const Eigen::ArrayXd>& x,
                       const Eigen::Ref<const Eigen::ArrayXd>& w) {
  const double n = static_cast<double>(x.size());
  const double a1 = -log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  dplog.resize(x.size());
  sqrt_neg_d2plog.resize(x.size());

  if (w.size() == 0) {
    for (int i = 0; i < x.size(); ++i) {
      if (n * x[i] < 1.0) {
        dplog[i] = a2 + 2.0 * a3 * x[i];
        sqrt_neg_d2plog[i] = a2 / 2.0;
        plog_sum += a1 + a2 * x[i] + a3 * x[i] * x[i];
      } else {
        dplog[i] = 1.0 / x[i];
        sqrt_neg_d2plog[i] = 1.0 / x[i];
        plog_sum += log(x[i]);
      }
    }
  } else {
    for (int i = 0; i < x.size(); ++i) {
      if (n * x[i] < w[i]) {
        dplog[i] = w[i] * (2.0 * n / w[i] - n * n * x[i] / (w[i] * w[i]));
        sqrt_neg_d2plog[i] = n / sqrt(w[i]);
        plog_sum += w[i] * (log(w[i] / n) - 1.5 + 2.0 * n *  x[i] / w[i] -
          0.5 * (n * n * x[i] * x[i]) / (w[i] * w[i]));
      } else {
        dplog[i] = w[i] / x[i];
        sqrt_neg_d2plog[i] = sqrt(w[i]) / x[i];
        plog_sum += w[i] * log(x[i]);
      }
    }
  }
}

Eigen::ArrayXd PSEUDO_LOG::plog(Eigen::VectorXd&& x) {
  const double n = static_cast<double>(x.size());
  const double a1 = -std::log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  for (int i = 0; i < x.size(); ++i) {
    if (n * x[i] < 1.0) {
      x[i] = a1 + a2 * x[i] + a3 * x[i] * x[i];
    } else {
      x[i] = log(x[i]);
    }
  }
  return x;
}

Eigen::ArrayXd PSEUDO_LOG::plog(Eigen::VectorXd&& x,
                                const Eigen::Ref<const Eigen::ArrayXd>& w) {
  const double n = static_cast<double>(x.size());
  const double a1 = -std::log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  Eigen::ArrayXd out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    if (n * x[i] < w[i]) {
      out[i] = (log(w[i] / n) - 1.5 + 2.0 * n *  x[i] / w[i] -
        0.5 * (n * n * x[i] * x[i]) / (w[i] * w[i]));
    } else {
      out[i] = log(x[i]);
    }
  }
  return out;
}

double PSEUDO_LOG::sum(const Eigen::Ref<const Eigen::VectorXd>& x,
                       const Eigen::Ref<const Eigen::ArrayXd>& w) {
  const double n = static_cast<double>(x.size());
  const double a1 = -log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  double out = 0;
  if (w.size() == 0) {
    for (int i = 0; i < x.size(); ++i) {
      out += n * x[i] < 1.0 ? a1 + a2 * x[i] + a3 * x[i] * x[i] : log(x[i]);
    }
  } else{
    for (int i = 0; i < x.size(); ++i) {
      out += n * x[i] < w[i] ?
      w[i] * (log(w[i] / n) - 1.5 + 2.0 * n *  x[i] / w[i] -
      0.5 * (n * n * x[i] * x[i]) / (w[i] * w[i])) :
      w[i] * log(x[i]);
    }
  }
  return out;
}
