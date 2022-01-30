#include "utils_lm.h"

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::VectorXd>& beta,
                     const Eigen::Ref<const Eigen::MatrixXd>& x,
                     const Eigen::Ref<const Eigen::VectorXd>& y) {
  return x.array().colwise() * (y - x * beta).array();
}

Eigen::VectorXd gradient_nlogLR_lm(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& x) {
  const Eigen::ArrayXd denominator =
    Eigen::VectorXd::Ones(g.rows()) + g * lambda;
  // gradient of nlogLR
  const Eigen::MatrixXd gradient =
    -(x.transpose() * (x.array().colwise() / denominator).matrix()) * lambda;
  return gradient;
}

minEL test_lm(const Eigen::Ref<const Eigen::VectorXd>& beta0,
              const Eigen::Ref<const Eigen::MatrixXd>& x,
              const Eigen::Ref<const Eigen::VectorXd>& y,
              const Eigen::Ref<const Eigen::MatrixXd>& lhs,
              const Eigen::Ref<const Eigen::VectorXd>& rhs,
              const double threshold,
              const int maxit,
              const double abstol) {
  /// initialization ///
  // orthogonal projection matrix
  const Eigen::MatrixXd P =
    Eigen::MatrixXd::Identity(lhs.cols(), lhs.cols()) -
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * lhs;
  // initial value (LS estimate) with the constraint imposed.
  Eigen::VectorXd beta =
    P * beta0 + lhs.transpose() * (lhs * lhs.transpose()).inverse() * rhs;
  // estimating function
  Eigen::MatrixXd g = g_lm(beta, x, y);
  // evaluation
  Eigen::VectorXd lambda = EL(g, threshold).lambda;
  // function value(-logLR)
  const int n = g.rows();
  double f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * lambda);

  /// minimization (projected gradient descent) ///
  double gamma = 1.0 / n;    // step size
  bool convergence = false;
  int iterations = 0;
  while (!convergence && iterations != maxit) {
    // update parameter
    Eigen::VectorXd beta_tmp =
      beta - gamma * P * gradient_nlogLR_lm(lambda, g, x);
    // update g
    Eigen::MatrixXd g_tmp = g_lm(beta_tmp, x, y);
    // update lambda
    EL eval(g_tmp, threshold);
    Eigen::VectorXd lambda_tmp = eval.lambda;
    if (!eval.convergence && iterations > 9) {
      return {beta, lambda, f1, iterations, convergence};
    }
    // update function value
    double f0 = f1;
    f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g_tmp * lambda_tmp);

    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < f1) {
      // reduce step size
      gamma /= 2;
      // propose new theta
      beta_tmp = beta - gamma * P * gradient_nlogLR_lm(lambda, g, x);;
      // propose new lambda
      g_tmp = g_lm(beta_tmp, x, y);
      lambda_tmp = EL(g_tmp, threshold).lambda;
      if (gamma < abstol) {
        return {beta, lambda, f0, iterations, convergence};
      }
      // propose new function value
      f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g_tmp * lambda_tmp);
    }

    // update
    beta = std::move(beta_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);
    ++iterations;

    // convergence check
    if (f0 - f1 < abstol) {
      convergence = true;
    }
  }
  return {beta, lambda, f1, iterations, convergence};
}
