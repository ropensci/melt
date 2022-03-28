#include "EL.h"

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd>&,
                              const Eigen::Ref<const Eigen::VectorXd>&)>
  EL::set_g_fcn(const std::string method)
{
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::MatrixXd>&,
             const Eigen::Ref<const Eigen::VectorXd>&)>>
               g_map{
               {{"mean", g_mean},
               {"lm", g_lm}}
             };
  return g_map[method];
}

std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::VectorXd>&,
                              const Eigen::Ref<const Eigen::MatrixXd>&,
                              const Eigen::Ref<const Eigen::MatrixXd>&)>
  EL::set_gr_fcn(const std::string method)
{
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> gr_map{
               {{"mean", gr_nloglr_mean},
                {"lm", gr_nloglr_lm}}
               };
  return gr_map[method];
}

void EL::set_el(const Eigen::Ref<const Eigen::MatrixXd>& g)
{
  // maximization
  while (!conv && iter != maxit && nllr <= th) {
    // pseudo log
    const PSEUDO_LOG pl(Eigen::VectorXd::Ones(n) + g * l);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * pl.sqrt_neg_d2plog;
    // propose new lambda by NR method with least square
    Eigen::VectorXd step = (J.transpose() * J).ldlt().solve(
      J.transpose() * (pl.dplog / pl.sqrt_neg_d2plog).matrix());
    // update function value
    nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step));
    // step halving to ensure increase in function value
    if (nllr < pl.plog_sum) {
      step /= 2;
      nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step));
    }
    // convergence check
    if (step.norm() < tol * l.norm() + tol) {
      conv = true;
    } else {
      ++iter;
    }
    // update lambda
    l += step;
  }
}

void EL::set_el(const Eigen::Ref<const Eigen::MatrixXd>& g,
                const Eigen::Ref<const Eigen::VectorXd>& w)
{
  // maximization
  while (!conv && iter != maxit && nllr <= th) {
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
    if (step.norm() < tol * l.norm() + tol) {
      conv = true;
    } else {
      ++iter;
    }
    // update lambda
    l += step;
  }
}

/* Constructor for EL class (evaluation)
 * Last updated: 03/21/21
 */
EL::EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
       const int maxit,
       const double tol,
       const double th)
  : l{Eigen::VectorXd::Zero(g.cols())},
    par{},
    maxit{maxit},
    tol{tol},
    th{th},
    n{static_cast<int>(g.rows())},
    g_fcn{},
    gr_fcn{}
{
  set_el(g);
}

EL::EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
       const Eigen::Ref<const Eigen::VectorXd>& w,
       const int maxit,
       const double tol,
       const double th)
  : l{Eigen::VectorXd::Zero(g.cols())},
    par{},
    maxit{maxit},
    tol{tol},
    th{th},
    n{static_cast<int>(g.rows())},
    g_fcn{},
    gr_fcn{}
{
  set_el(g, w);
}

EL::EL(const std::string method,
       const Eigen::Ref<const Eigen::VectorXd>& par0,
       const Eigen::Ref<const Eigen::MatrixXd>& x,
       const int maxit,
       const double tol,
       const double th)
  : l{Eigen::VectorXd::Zero(x.cols())},
    par{par0},
    maxit{maxit},
    tol{tol},
    th{th},
    n{static_cast<int>(x.rows())},
    g_fcn{set_g_fcn(method)},
    gr_fcn{set_gr_fcn(method)}
{
  set_el(g_fcn(x, par));
}

EL::EL(const std::string method,
       const Eigen::Ref<const Eigen::VectorXd>& par0,
       const Eigen::Ref<const Eigen::MatrixXd>& x,
       const Eigen::Ref<const Eigen::VectorXd>& w,
       const int maxit,
       const double tol,
       const double th)
  : l{Eigen::VectorXd::Zero(x.cols())},
    par{par0},
    maxit{maxit},
    tol{tol},
    th{th},
    n{static_cast<int>(x.rows())},
    g_fcn{set_g_fcn(method)},
    gr_fcn{set_gr_fcn(method)}
{
  set_el(g_fcn(x, par), w);
}

/* Constructor for EL class (minimization)
 * Last updated: 03/21/21
 */
EL::EL(const std::string method,
       const Eigen::Ref<const Eigen::VectorXd>& par0,
       const Eigen::Ref<const Eigen::MatrixXd>& x,
       const Eigen::Ref<const Eigen::MatrixXd>& lhs,
       const Eigen::Ref<const Eigen::VectorXd>& rhs,
       const int maxit,
       const double tol,
       const double th)
  : maxit{maxit},
    tol{tol},
    th{th},
    n{static_cast<int>(x.rows())},
    g_fcn{set_g_fcn(method)},
    gr_fcn{set_gr_fcn(method)}
{
  /// initialization ///
  // orthogonal projection matrix
  const Eigen::MatrixXd proj =
    Eigen::MatrixXd::Identity(lhs.cols(), lhs.cols()) -
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * lhs;
  // parameter (constraint imposed)
  par = proj * par0 + lhs.transpose() * (lhs * lhs.transpose()).inverse() * rhs;
  // estimating function
  Eigen::MatrixXd g = g_fcn(x, par);
  // lambda
  l = EL(g, maxit, tol, th).l;
  // function value(-logLR)
  nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * l);
  // function norm
  const double norm0 = (proj * gr_fcn(l, g, x)).norm();

  /// minimization (projected gradient descent) ///
  double gamma = 1.0;
  while (!conv && iter != maxit) {
    // update parameter
    Eigen::VectorXd par_tmp = par - gamma * proj * gr_fcn(l, g, x);
    // update estimating function
    Eigen::MatrixXd g_tmp = g_fcn(x, par_tmp);
    // update lambda
    Eigen::VectorXd l_tmp = EL(g_tmp, maxit, tol, th).l;
    // update function value
    double f0 = nllr;
    nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g_tmp * l_tmp);
    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < nllr) {
      // reduce step size
      gamma /= 2.0;
      // propose new parameter
      par_tmp = par - gamma * proj * gr_fcn(l, g, x);
      // propose new lambda
      g_tmp = g_fcn(x, par_tmp);
      l_tmp = EL(g_tmp, maxit, tol, th).l;
      if (gamma < tol) {
        nllr = f0;
        break;
      }
      // propose new function value
      nllr = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g_tmp * l_tmp);
    }
    // update
    par = std::move(par_tmp);
    l = std::move(l_tmp);
    g = std::move(g_tmp);
    // convergence check
    if ((proj * gr_fcn(l, g, x)).norm() < tol * norm0 + tol) {
      conv = true;
    } else {
      ++iter;
    }
  }
}

/* log probability for weighted EL class
 * Last updated: 03/09/21
 *
 */
Eigen::ArrayXd EL::log_prob(const Eigen::Ref<const Eigen::MatrixXd>& x,
                            const Eigen::Ref<const Eigen::ArrayXd>& w) const {
  const Eigen::MatrixXd g = g_fcn(x, par);
  return  w.log() - PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g * l);
}

/* log weighted probability for weighted EL class
 * Last updated: 03/09/21
 *
 */
Eigen::ArrayXd EL::log_wprob(const Eigen::Ref<const Eigen::MatrixXd>& x,
                             const Eigen::Ref<const Eigen::ArrayXd>& w) const {
  const Eigen::MatrixXd g = g_fcn(x, par);
  return  -w * PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g * l);
}


// Eigen::ArrayXd EL::logp(const Eigen::Ref<const Eigen::MatrixXd>& x,
//                         const Eigen::Ref<const Eigen::ArrayXd>& w) {
//   const Eigen::MatrixXd g = g_fcn(x, par);
//   return  w.log() - PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g * l);
// }







/* Constructor for PSEUDO_LOG class
 * Last updated: 03/19/21
 *
 */
PSEUDO_LOG::PSEUDO_LOG(Eigen::VectorXd&& x) {
  const double n = static_cast<double>(x.size());
  const double a1 = -log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;

  dplog.resize(x.size());
  sqrt_neg_d2plog.resize(x.size());

  for (unsigned int i = 0; i < x.size(); ++i) {
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
}

/* Constructor for PSEUDO_LOG class (weighted)
 * Last updated: 03/19/21
 *
 */
PSEUDO_LOG::PSEUDO_LOG(Eigen::VectorXd&& x,
                       const Eigen::Ref<const Eigen::VectorXd>& w) {
  const double n = static_cast<double>(x.size());
  const double a1 = -log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;

  dplog.resize(x.size());
  sqrt_neg_d2plog.resize(x.size());

  for (unsigned int i = 0; i < x.size(); ++i) {
    if (n * x[i] < 1.0) {
      dplog[i] = w[i] * (a2 + 2.0 * a3 * x[i]);
      sqrt_neg_d2plog[i] = sqrt(w[i]) * a2 / 2.0;
      plog_sum += w[i] * (a1 + a2 * x[i] + a3 * x[i] * x[i]);
    } else {
      dplog[i] = w[i] / x[i];
      sqrt_neg_d2plog[i] = sqrt(w[i]) / x[i];
      plog_sum += w[i] * log(x[i]);
    }
  }
}

/* pseudo log function
 * Last updated: 03/19/21
 *
 */
Eigen::ArrayXd PSEUDO_LOG::plog(Eigen::VectorXd&& x) {
  const double n = static_cast<double>(x.size());
  const double a1 = -log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  for (unsigned int i = 0; i < x.size(); ++i) {
    if (n * x[i] < 1.0) {
      x[i] = a1 + a2 * x[i] + a3 * x[i] * x[i];
    } else {
      x[i] = log(x[i]);
    }
  }
  return x;
}

double PSEUDO_LOG::sum(Eigen::VectorXd&& x) {
  const double n = static_cast<double>(x.size());
  const double a1 = -log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  double out = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    out += n * x[i] < 1.0 ? a1 + a2 * x[i] + a3 * x[i] * x[i] : log(x[i]);
  }
  return out;
}

/* Summation function for PSEUDO_LOG class (weighted)
 * Last updated: 03/19/21
 *
 */
double PSEUDO_LOG::sum(Eigen::VectorXd&& x,
                       const Eigen::Ref<const Eigen::VectorXd>& w) {
  const double n = static_cast<double>(x.size());
  const double a1 = -log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = -0.5 * n * n;
  double out = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    out += n * x[i] < 1.0 ?
    w[i] * (a1 + a2 * x[i] + a3 * x[i] * x[i]) :
    w[i] * log(x[i]);
  }
  return out;
}
