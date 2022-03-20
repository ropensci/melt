#include "EL.h"

/* Constructor for EL class (evaluation)
 * Last updated: 03/19/21
 */
EL::EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
       const int maxit,
       const double tol,
       const double threshold)
  : n{static_cast<int>(g.rows())},
    l{Eigen::VectorXd::Zero(g.cols())}
{
  // maximization
  while (!convergence && iterations != maxit && nlogLR <= threshold) {
    // pseudo log
    const PSEUDO_LOG pl(Eigen::VectorXd::Ones(n) + g * l);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * pl.sqrt_neg_d2plog;
    // propose new lambda by NR method with least square
    Eigen::VectorXd step = (J.transpose() * J).ldlt().solve(
      J.transpose() * (pl.dplog / pl.sqrt_neg_d2plog).matrix());
    // update function value
    nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step));
    // step halving to ensure increase in function value
    if (nlogLR < pl.plog_sum) {
      step /= 2;
      nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step));
    }
    // convergence check
    if (step.norm() < tol * l.norm() + tol) {
      convergence = true;
    } else {
      ++iterations;
    }
    // update lambda
    l += step;
  }
}

/* Constructor for weighted EL class (evaluation)
 * Last updated: 03/19/21
 */
EL::EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
       const Eigen::Ref<const Eigen::VectorXd>& w,
       const int maxit,
       const double tol,
       const double threshold)
  : n{static_cast<int>(g.rows())},
    l{Eigen::VectorXd::Zero(g.cols())}
{
  // maximization
  while (!convergence && iterations != maxit && nlogLR <= threshold) {
    // pseudo log
    const PSEUDO_LOG pl(Eigen::VectorXd::Ones(n) + g * l, w);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * pl.sqrt_neg_d2plog;
    // propose new lambda by NR method with least square
    Eigen::VectorXd step = (J.transpose() * J).ldlt().solve(
      J.transpose() * (pl.dplog / pl.sqrt_neg_d2plog).matrix());
    // update function value
    nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step), w);
    // step halving to ensure increase in function value
    if (nlogLR < pl.plog_sum) {
      step /= 2;
      nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (l + step), w);
    }
    // convergence check
    if (step.norm() < tol * l.norm() + tol) {
      convergence = true;
    } else {
      ++iterations;
    }
    // update lambda
    l += step;
  }
}

/* Constructor for EL2 class (evaluation)
 * Last updated: 03/09/21
 *
 * abstol for gamma should be reconsidered.
 * Perhaps, use another optim parameter such as step size tolerance.
 */
EL2::EL2(const std::string method,
         const Eigen::Ref<const Eigen::VectorXd>& par0,
         const Eigen::Ref<const Eigen::MatrixXd>& x,
         const int maxit,
         const double abstol,
         const double threshold)
  : type{method}, par{par0}, n{static_cast<int>(x.rows())}
{
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> funcMap{
               {{"mean", g_mean},
               {"lm", g_lm}}
             };
  Eigen::MatrixXd g = funcMap[type](par, x);
  // maximization
  lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  // lambda = Eigen::VectorXd::Zero(par0.size());
  while (!convergence && iterations != maxit && nlogLR <= threshold) {
    // plog class
    PSEUDO_LOG pl(Eigen::VectorXd::Ones(n) + g * lambda);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * pl.sqrt_neg_d2plog;
    // propose new lambda by NR method with least square
    Eigen::VectorXd step = (J.transpose() * J).ldlt().solve(
      J.transpose() * (pl.dplog / pl.sqrt_neg_d2plog).matrix());
    // update function value
    nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (lambda + step));
    // step halving to ensure increase in function value
    double gamma = 1.0;
    while (nlogLR < pl.plog_sum) {
      gamma /= 2;
      if (gamma < abstol) {
        break;
      }
      nlogLR =
        PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (lambda + gamma * step));
    }
    /* If the step halving is not successful (possibly due to the convex
     * hull constraint), terminate the maximization with the current values
     * without further updates.
     */
    if (gamma < abstol) {
      nlogLR = pl.plog_sum;
      break;
    }
    // Otherwise, update lambda and check for convergence
    lambda += gamma * step;
    if (nlogLR - pl.plog_sum < abstol) {
      convergence = true;
    } else {
      ++iterations;
    }
  }
}

/* Constructor for weighted EL2 class (evaluation)
 * Last updated: 03/19/21
 *
 * abstol for gamma should be reconsidered.
 * Perhaps, use another optim parameter such as step size tolerance.
 *
 * l2 norm of the gradient of nlogLR is
 * (g.array().colwise() * pl.dplog).colwise().sum().matrix().norm(), where
 * the pl is computed with the current value of lambda. Consider moving
 * toward the real gradient descent than evaluating nlogLR directly. It could be
 * slower.
 */
EL2::EL2(const std::string method,
         const Eigen::Ref<const Eigen::VectorXd>& par0,
         const Eigen::Ref<const Eigen::MatrixXd>& x,
         const Eigen::Ref<const Eigen::VectorXd>& w,
         const int maxit,
         const double abstol,
         const double threshold)
  : type{method}, par{par0}, n{static_cast<int>(x.rows())}
{
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> funcMap{
               {{"mean", g_mean},
               {"lm", g_lm}}
             };
  Eigen::MatrixXd g = funcMap[type](par, x);
  // maximization
  lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  while (!convergence && iterations != maxit && nlogLR <= threshold) {
    // plog class
    PSEUDO_LOG pl(Eigen::VectorXd::Ones(n) + g * lambda, w);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * pl.sqrt_neg_d2plog;
    // propose new lambda by NR method with least square
    Eigen::VectorXd step = (J.transpose() * J).ldlt().solve(
      J.transpose() * (pl.dplog / pl.sqrt_neg_d2plog).matrix());
    // update function value
    nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (lambda + step), w);
    // step halving to ensure increase in function value
    double gamma = 1.0;
    while (nlogLR < pl.plog_sum) {
      gamma /= 2;
      if (gamma < abstol) {
        break;
      }
      nlogLR =
        PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (lambda + gamma * step),
                        w);
    }
    /* If the step halving is not successful (possibly due to the convex
     * hull constraint), terminate the maximization with the current values
     * without further updates.
     */
    if (gamma < abstol) {
      nlogLR = pl.plog_sum;
      break;
    }
    // Otherwise, update lambda and check for convergence
    lambda += gamma * step;
    if (nlogLR - pl.plog_sum < abstol) {
      convergence = true;
    } else {
      ++iterations;
    }
  }
}

/* Constructor for EL2 class (minimization)
 * Last updated: 03/09/21
 */
EL2::EL2(const std::string method,
         const Eigen::Ref<const Eigen::VectorXd>& par0,
         const Eigen::Ref<const Eigen::MatrixXd>& x,
         const Eigen::Ref<const Eigen::MatrixXd>& lhs,
         const Eigen::Ref<const Eigen::VectorXd>& rhs,
         const int maxit,
         const double abstol,
         const double threshold)
  : type{method}, par{par0}, n{static_cast<int>(x.rows())}
{
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> funcMap{
               {{"mean", g_mean},
               {"lm", g_lm}}
             };
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> gradMap{
               {{"mean", gr_nlogLR_lm},
               {"lm", gr_nlogLR_lm}}
             };
  /// initialization ///
  // orthogonal projection matrix
  const Eigen::MatrixXd P =
    Eigen::MatrixXd::Identity(lhs.cols(), lhs.cols()) -
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * lhs;
  // initial value (LS estimate) with the constraint imposed.
  par = P * par0 + lhs.transpose() * (lhs * lhs.transpose()).inverse() * rhs;
  // estimating function
  Eigen::MatrixXd g = funcMap[type](par, x);
  // evaluation
  lambda = EL(g, maxit, abstol, threshold).l;
  // function value(-logLR)
  nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * lambda);

  /// minimization (projected gradient descent) ///
  double gamma = 1.0 / n;    // step size
  convergence = false;
  iterations = 1;
  while (!convergence && iterations != maxit) {
    // update parameter
    Eigen::VectorXd par_tmp =
      par - gamma * P * gradMap[type](lambda, g, x);
    // update g
    Eigen::MatrixXd g_tmp = funcMap[type](par_tmp, x);
    // update lambda
    EL eval(g_tmp, maxit, abstol, threshold);
    Eigen::VectorXd lambda_tmp = eval.l;
    if (!eval.convergence && iterations > 9) {
      break;
    }

    // update function value
    double f0 = nlogLR;
    nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g_tmp * lambda_tmp);

    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < nlogLR) {
      // reduce step size
      gamma /= 2.0;
      // propose new par
      par_tmp = par - gamma * P * gradMap[type](lambda, g, x);
      // propose new lambda
      g_tmp = funcMap[type](par_tmp, x);
      lambda_tmp = EL(g_tmp, maxit, abstol, threshold).l;
      if (gamma < abstol) {
        nlogLR = f0;
        break;
      }
      // propose new function value
      nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g_tmp * lambda_tmp);
    }

    // update
    par = std::move(par_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);
    ++iterations;

    // convergence check
    if (f0 - nlogLR < abstol) {
      convergence = true;
    }
  }
}

/* log probability for weighted EL2 class
 * Last updated: 03/09/21
 *
 */
Eigen::ArrayXd EL2::log_prob(const Eigen::Ref<const Eigen::MatrixXd>& x,
                             const Eigen::Ref<const Eigen::ArrayXd>& w) {
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> funcMap{
               {{"mean", g_mean},
               {"lm", g_lm}}
             };
  Eigen::MatrixXd g = funcMap[type](par, x);
  return  w.log() -
    PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g * lambda);
}

/* log weighted probability for weighted EL2 class
 * Last updated: 03/09/21
 *
 */
Eigen::ArrayXd EL2::log_wprob(const Eigen::Ref<const Eigen::MatrixXd>& x,
                              const Eigen::Ref<const Eigen::ArrayXd>& w) {
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> funcMap{
               {{"mean", g_mean},
               {"lm", g_lm}}
             };
  Eigen::MatrixXd g = funcMap[type](par, x);
  return  -w * PSEUDO_LOG::plog(Eigen::VectorXd::Ones(n) + g * lambda);
}





/* Constructor for PSEUDO_LOG class
 * Last updated: 03/04/21
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

double th_nlogLR(const int p, const Rcpp::Nullable<double> threshold) {
  return (threshold.isNull())? 20.0 * p : Rcpp::as<double>(threshold);
};

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::VectorXd>& par,
                       const Eigen::Ref<const Eigen::MatrixXd>& x) {
  return x.rowwise() - par.transpose();
}

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::VectorXd>& par,
                     const Eigen::Ref<const Eigen::MatrixXd>& data) {
  // const Eigen::VectorXd y = data.col(0);
  // const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  // return x.array().colwise() * (y - x * beta).array();
  return data.rightCols(data.cols() - 1).array().colwise() *
    (data.col(0) - data.rightCols(data.cols() - 1) * par).array();
}

Eigen::VectorXd gr_nlogLR_lm(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data) {
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd denominator =
    Eigen::VectorXd::Ones(g.rows()) + g * lambda;
  // const Eigen::MatrixXd gradient =
  //   -(x.transpose() * (x.array().colwise() / denominator).matrix()) * lambda;
  return
    -(x.transpose() * (x.array().colwise() / denominator).matrix()) * lambda;
}
