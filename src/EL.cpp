#include "EL.h"
#include "helpers.h"
#include <RcppEigen.h>
#include <cfloat>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <utility>

/* EL class (evaluation)
 * Last updated: 07/07/22
 */
EL::EL(const std::string method,
       const Eigen::Ref<const Eigen::VectorXd> &par0,
       const Eigen::Ref<const Eigen::MatrixXd> &x,
       const int maxit_l,
       const double tol_l,
       const double th,
       const Eigen::Ref<const Eigen::ArrayXd> &wt)
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

EL::EL(const Eigen::Ref<const Eigen::MatrixXd> &g,
       const int maxit_l,
       const double tol_l,
       const double th,
       const Eigen::Ref<const Eigen::ArrayXd> &wt)
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

std::function<Eigen::VectorXd(const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::ArrayXd> &)>
EL::set_mele_fn(const std::string method)
{
  std::map<std::string, std::function<Eigen::VectorXd(
                            const Eigen::Ref<const Eigen::MatrixXd> &,
                            const Eigen::Ref<const Eigen::ArrayXd> &)>>
      mele_map{{{"mean", mele_mean},
                {"sd", mele_sd},
                {"lm", mele_lm}}};
  return mele_map[method];
}

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &)>
EL::set_g_fn(const std::string method)
{
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
             {"quasipoisson_identity", g_qpoi_identity}}};
  return g_map[method];
}

void EL::set_el(const Eigen::Ref<const Eigen::MatrixXd> &g,
                const Eigen::Ref<const Eigen::ArrayXd> &w)
{
  // maximization
  while (!conv && iter != maxit_l && nllr <= th)
  {
    // pseudo log
    const PseudoLog pl(Eigen::VectorXd::Ones(n) + g * l, w);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * pl.sqrt_neg_d2plog;
    // propose new lambda by NR method with least square
    Eigen::VectorXd step = J.householderQr().solve(
        (pl.dplog / pl.sqrt_neg_d2plog).matrix());
    // update function value
    nllr = PseudoLog::sum(Eigen::VectorXd::Ones(n) + g * (l + step), w);
    // step halving to ensure increase in function value
    if (nllr < pl.plog_sum)
    {
      step /= 2;
      nllr = PseudoLog::sum(Eigen::VectorXd::Ones(n) + g * (l + step), w);
    }
    // convergence check
    if (step.norm() < tol_l * l.norm() + tol_l * tol_l)
    {
      conv = true;
    }
    ++iter;
    // update lambda
    l += step;
  }
}

Eigen::ArrayXd EL::logp_g(const Eigen::Ref<const Eigen::MatrixXd> &g) const
{
  if (w.size() == 0)
  {
    return -std::log(n) - PseudoLog::plog(Eigen::VectorXd::Ones(n) + g * l);
  }
  else
  {
    return log(w) - std::log(n) -
           PseudoLog::plog(Eigen::VectorXd::Ones(n) + g * l, w);
  }
}

Eigen::ArrayXd EL::logp(const Eigen::Ref<const Eigen::MatrixXd> &x) const
{
  if (w.size() == 0)
  {
    return -std::log(n) -
           PseudoLog::plog(Eigen::VectorXd::Ones(n) + g_fn(x, par) * l);
  }
  else
  {
    return log(w) - std::log(n) -
           PseudoLog::plog(Eigen::VectorXd::Ones(n) + g_fn(x, par) * l, w);
  }
}

double EL::loglik() const
{
  if (w.size() == 0)
  {
    return -nllr - n * std::log(n);
  }
  else
  {
    return -nllr - (w * (std::log(n) - log(w))).sum();
  }
}


/* CEL class (minimization)
 * Last updated: 08/23/22
 */
CEL::CEL(const std::string method,
         const Eigen::Ref<const Eigen::VectorXd> &par0,
         const Eigen::Ref<const Eigen::MatrixXd> &x,
         const Eigen::Ref<const Eigen::MatrixXd> &lhs,
         const Eigen::Ref<const Eigen::VectorXd> &rhs,
         const int maxit,
         const int maxit_l,
         const double tol,
         const double tol_l,
         const double step,
         const double th,
         const Eigen::Ref<const Eigen::ArrayXd> &wt)
    : gamma{step},
      n{static_cast<int>(x.rows())},
      weighted{wt.size() != 0},
      g_fn{CEL::set_g_fn(method)},
      gr_fn{CEL::set_gr_fn(method)}
{
  /// initialization ///
  // parameter (constraint imposed)
  par = proj(lhs, par0) +
        lhs.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
  // estimating function
  Eigen::MatrixXd g = g_fn(x, par);
  // lambda
  l = EL(g, maxit_l, tol_l, th, wt).l;
  // function value (-logLR)
  nllr = PseudoLog::sum(Eigen::VectorXd::Ones(n) + g * l, wt);
  // function norm
  // const double norm0 = proj(lhs, gr_fn(l, g, x, par, wt, weighted)).norm();

  /// minimization (projected gradient descent) ///
  while (!conv && iter != maxit && nllr <= th)
  {
    // update parameter
    Eigen::VectorXd par_tmp =
        par - gamma * proj(lhs, gr_fn(l, g, x, par, wt, weighted));
    // update estimating function
    Eigen::MatrixXd g_tmp = g_fn(x, par_tmp);
    // update lambda
    Eigen::VectorXd l_tmp = EL(g_tmp, maxit_l, tol_l, th, wt).l;
    // update function value
    const double f0 = nllr;
    nllr = PseudoLog::sum(Eigen::VectorXd::Ones(n) + g_tmp * l_tmp, wt);
    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < nllr || nllr < 0)
    {
      // reduce step size
      gamma /= 2;
      if (gamma < DBL_EPSILON)
      {
        break;
      }
      // propose new parameter
      par_tmp = par - gamma * proj(lhs, gr_fn(l, g, x, par, wt, weighted));
      // propose new lambda
      g_tmp = g_fn(x, par_tmp);
      l_tmp = EL(g_tmp, maxit_l, tol_l, th, wt).l;
      // propose new function value
      nllr = PseudoLog::sum(Eigen::VectorXd::Ones(n) + g_tmp * l_tmp, wt);
    }
    if (nllr >= f0)
    {
      nllr = f0;
      ++iter;
      break;
    }
    else if (std::isnan(nllr))
    {
      // return initial values
      par = proj(lhs, par0) +
            lhs.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
      Eigen::MatrixXd g = g_fn(x, par);
      l = EL(g, maxit_l, tol_l, th, wt).l;
      nllr = PseudoLog::sum(Eigen::VectorXd::Ones(n) + g * l, wt);
    }
    // update
    const double s = (par - par_tmp).norm();
    const double d = par.norm();
    par = std::move(par_tmp);
    l = std::move(l_tmp);
    g = std::move(g_tmp);
    // convergence check
    if (proj(lhs, gr_fn(l, g, x, par, wt, weighted)).norm() < tol ||
        s < tol * d + tol * tol)
    {
      conv = true;
    }
    ++iter;
  }
}

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &)>
CEL::set_g_fn(const std::string method)
{
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
             {"quasipoisson_identity", g_qpoi_identity}}};
  return g_map[method];
}

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::VectorXd> &,
                              const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::MatrixXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &,
                              const Eigen::Ref<const Eigen::ArrayXd> &,
                              const bool)>
CEL::set_gr_fn(const std::string method)
{
  std::map<std::string, std::function<Eigen::MatrixXd(
                            const Eigen::Ref<const Eigen::VectorXd> &,
                            const Eigen::Ref<const Eigen::MatrixXd> &,
                            const Eigen::Ref<const Eigen::MatrixXd> &,
                            const Eigen::Ref<const Eigen::VectorXd> &,
                            const Eigen::Ref<const Eigen::ArrayXd> &,
                            const bool)>>
      gr_map{
          {{"mean", gr_nloglr_mean},
           {"sd", gr_nloglr_sd},
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
           {"quasipoisson_log", gr_nloglr_qpoi_log},
           {"quasipoisson_identity", gr_nloglr_qpoi_identity}}};
  return gr_map[method];
}

Eigen::ArrayXd CEL::logp(const Eigen::Ref<const Eigen::MatrixXd> &x,
                         const Eigen::Ref<const Eigen::ArrayXd> &wt) const
{
  if (weighted)
  {
    return log(wt) - std::log(n) -
           PseudoLog::plog(Eigen::VectorXd::Ones(n) + g_fn(x, par) * l, wt);
  }
  else
  {
    return -std::log(n) -
           PseudoLog::plog(Eigen::VectorXd::Ones(n) + g_fn(x, par) * l);
  }
}

double CEL::loglik(const Eigen::Ref<const Eigen::ArrayXd> &wt) const
{
  if (weighted)
  {
    return -nllr - (wt * (std::log(n) - log(wt))).sum();
  }
  else
  {
    return -nllr - n * std::log(n);
  }
}


PseudoLog::PseudoLog(const Eigen::Ref<const Eigen::ArrayXd> &x,
                     const Eigen::Ref<const Eigen::ArrayXd> &w)
{
  const double n = static_cast<double>(x.size());
  const double a1 = -std::log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  dplog.resize(x.size());
  sqrt_neg_d2plog.resize(x.size());

  if (w.size() == 0)
  {
    for (int i = 0; i < x.size(); ++i)
    {
      if (n * x[i] < 1.0)
      {
        dplog[i] = a2 + 2.0 * a3 * x[i];
        sqrt_neg_d2plog[i] = a2 / 2.0;
        plog_sum += a1 + a2 * x[i] + a3 * x[i] * x[i];
      }
      else
      {
        dplog[i] = 1.0 / x[i];
        sqrt_neg_d2plog[i] = 1.0 / x[i];
        plog_sum += std::log(x[i]);
      }
    }
  }
  else
  {
    for (int i = 0; i < x.size(); ++i)
    {
      if (n * x[i] < w[i])
      {
        dplog[i] = w[i] * (a2 / w[i] - n * n * x[i] / (w[i] * w[i]));
        sqrt_neg_d2plog[i] = n / sqrt(w[i]);
        plog_sum += w[i] * (log(w[i] / n) - 1.5 + a2 * x[i] / w[i] +
                            a3 * (x[i] * x[i]) / (w[i] * w[i]));
      }
      else
      {
        dplog[i] = w[i] / x[i];
        sqrt_neg_d2plog[i] = sqrt(w[i]) / x[i];
        plog_sum += w[i] * std::log(x[i]);
      }
    }
  }
}

Eigen::ArrayXd PseudoLog::plog(Eigen::VectorXd &&x)
{
  const double n = static_cast<double>(x.size());
  const double a1 = -std::log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  for (int i = 0; i < x.size(); ++i)
  {
    if (n * x[i] < 1.0)
    {
      x[i] = a1 + a2 * x[i] + a3 * x[i] * x[i];
    }
    else
    {
      x[i] = std::log(x[i]);
    }
  }
  return x;
}

Eigen::ArrayXd PseudoLog::plog(Eigen::VectorXd &&x,
                               const Eigen::Ref<const Eigen::ArrayXd> &w)
{
  const double n = static_cast<double>(x.size());
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  Eigen::ArrayXd out(x.size());
  for (int i = 0; i < x.size(); ++i)
  {
    if (n * x[i] < w[i])
    {
      out[i] = log(w[i] / n) - 1.5 + a2 * x[i] / w[i] +
               a3 * (x[i] * x[i]) / (w[i] * w[i]);
    }
    else
    {
      out[i] = log(x[i]);
    }
  }
  return out;
}

double PseudoLog::sum(const Eigen::Ref<const Eigen::VectorXd> &x,
                      const Eigen::Ref<const Eigen::ArrayXd> &w)
{
  const double n = static_cast<double>(x.size());
  const double a1 = -std::log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  double out = 0;
  if (w.size() == 0)
  {
    for (int i = 0; i < x.size(); ++i)
    {
      out +=
          n * x[i] < 1.0 ? a1 + a2 * x[i] + a3 * x[i] * x[i] : std::log(x[i]);
    }
  }
  else
  {
    for (int i = 0; i < x.size(); ++i)
    {
      out +=
          n * x[i] < w[i] ? w[i] *
                                (std::log(w[i] / n) - 1.5 + a2 * x[i] / w[i] +
                                 a3 * (x[i] * x[i]) / (w[i] * w[i]))
                          : w[i] * std::log(x[i]);
    }
  }
  return out;
}
