#include "EL.h"

EL::EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
       const double threshold,
       const int maxit,
       const double abstol) {
  // maximization
  lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  iterations = 1;
  convergence = false;
  while (!convergence && iterations != maxit) {
    // plog class
    PSEUDO_LOG log_tmp(Eigen::VectorXd::Ones(g.rows()) + g * lambda);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * log_tmp.sqrt_neg_d2plog;
    // prpose new lambda by NR method with least square
    Eigen::VectorXd step =
      (J.transpose() * J).ldlt().solve(
          J.transpose() * (log_tmp.dplog / log_tmp.sqrt_neg_d2plog).matrix());
    // update function value
    nlogLR =
      PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
    // step halving to ensure validity
    while (nlogLR < log_tmp.plog_sum) {
      step /= 2;
      nlogLR =
        PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
    }
    // update lambda
    lambda += step;

    // check convex hull constraint(stop if larger than threshold)
    if (nlogLR > threshold) {
      break;
    }

    // convergence check
    if (nlogLR - log_tmp.plog_sum < abstol) {
      convergence = true;
    } else {
      ++iterations;
    }
  }
}

EL2::EL2(const Eigen::Ref<const Eigen::VectorXd>& par0,
         const Eigen::Ref<const Eigen::MatrixXd>& x,
         const std::string type,
         const double threshold,
         const int maxit,
         const double abstol) {
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> funcMap{
               {{"mean", g_mean},
               {"lm", g_lm2}}
             };

  Eigen::MatrixXd g = funcMap[type](par0, x);
  // maximization
  lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  iterations = 1;
  convergence = false;
  while (!convergence && iterations != maxit) {
    // plog class
    PSEUDO_LOG log_tmp(Eigen::VectorXd::Ones(g.rows()) + g * lambda);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * log_tmp.sqrt_neg_d2plog;
    // prpose new lambda by NR method with least square
    Eigen::VectorXd step =
      (J.transpose() * J).ldlt().solve(
          J.transpose() * (log_tmp.dplog / log_tmp.sqrt_neg_d2plog).matrix());
    // update function value
    nlogLR =
      PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
    // step halving to ensure validity
    while (nlogLR < log_tmp.plog_sum) {
      step /= 2;
      nlogLR =
        PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
    }
    // update lambda
    lambda += step;

    // check convex hull constraint(stop if larger than threshold)
    if (nlogLR > threshold) {
      break;
    }

    // convergence check
    if (nlogLR - log_tmp.plog_sum < abstol) {
      convergence = true;
    } else {
      ++iterations;
    }
  }
}

EL2::EL2(const Eigen::Ref<const Eigen::VectorXd>& par0,
         const Eigen::Ref<const Eigen::MatrixXd>& x,
         const std::string type,
         const Eigen::Ref<const Eigen::MatrixXd>& lhs,
         const Eigen::Ref<const Eigen::VectorXd>& rhs,
         const double threshold,
         const int maxit,
         const double abstol) {
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> funcMap{
               {{"mean", g_mean},
               {"lm", g_lm2}}
             };
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> gradMap{
               {{"mean", gradient_nlogLR_lm2},
               {"lm", gradient_nlogLR_lm2}}
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
  lambda = EL(g, threshold).lambda;
  // function value(-logLR)
  const int n = g.rows();
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
    EL eval(g_tmp, threshold);
    Eigen::VectorXd lambda_tmp = eval.lambda;
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
      lambda_tmp = EL(g_tmp, threshold).lambda;
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

PSEUDO_LOG::PSEUDO_LOG(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -std::log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;

  dplog.resize(x.size());
  sqrt_neg_d2plog.resize(x.size());
  plog_sum = 0;

  for (unsigned int i = 0; i < x.size(); ++i) {
    if (x[i] < a0) {
      dplog[i] = a2 + 2 * a3 * x[i];
      sqrt_neg_d2plog[i] = a2 / 2;
      plog_sum += a1 + a2 * x[i] + a3 * x[i] * x[i];
    } else {
      dplog[i] = 1.0 / x[i];
      sqrt_neg_d2plog[i] = 1.0 / x[i];
      plog_sum += std::log(x[i]);
    }
  }
}

double PSEUDO_LOG::sum(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -std::log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;
  double out = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    out += x[i] < a0 ? a1 + a2 * x[i] + a3 * x[i] * x[i] : std::log(x[i]);
  }
  return out;
}

Eigen::ArrayXd PSEUDO_LOG::dp(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = 2.0 * n;
  static const double a2 = -1.0 * n * n;
  // Eigen::ArrayXd out(n);
  for (unsigned int i = 0; i < x.size(); ++i) {
    if (x[i] < a0) {
      x[i] = a1 + a2 * x[i];
    } else {
      x[i] = 1.0 / x[i];
    }
  }
  return x;
}

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::VectorXd>& par,
                       const Eigen::Ref<const Eigen::MatrixXd>& x) {
  return x.rowwise() - par.transpose();
}

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::VectorXd>& beta,
                     const Eigen::Ref<const Eigen::MatrixXd>& x,
                     const Eigen::Ref<const Eigen::VectorXd>& y) {
  return x.array().colwise() * (y - x * beta).array();
}

Eigen::MatrixXd g_lm2(const Eigen::Ref<const Eigen::VectorXd>& par,
                      const Eigen::Ref<const Eigen::MatrixXd>& data) {
  // const Eigen::VectorXd y = data.col(0);
  // const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  // return x.array().colwise() * (y - x * beta).array();
  return data.rightCols(data.cols() - 1).array().colwise() *
    (data.col(0) - data.rightCols(data.cols() - 1) * par).array();
}

Eigen::VectorXd gradient_nlogLR_lm2(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data) {
  const Eigen::MatrixXd x = data.rightCols(data.cols() - 1);
  const Eigen::ArrayXd denominator =
    Eigen::VectorXd::Ones(g.rows()) + g * lambda;
  // gradient of nlogLR
  const Eigen::MatrixXd gradient =
    -(x.transpose() * (x.array().colwise() / denominator).matrix()) * lambda;
    return gradient;
}
