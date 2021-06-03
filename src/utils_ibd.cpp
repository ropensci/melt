#include "utils_ibd.h"

arma::mat g_ibd(const arma::vec& theta,
                const arma::mat& x,
                const arma::mat& c) {
  return x - c.each_row() % theta.t();
}


arma::mat cov_ibd(const arma::mat& x,
                  const arma::mat& c,
                  const bool adjust) {
  // number of blocks
  const int n = x.n_rows;
  // estimator(global minimizer)
  arma::vec theta_hat = n * arma::trans(arma::mean(x, 0) / arma::sum(c, 0));
  // estimating function
  arma::mat g = x - c.each_row() % theta_hat.t();
  // covariance estimate(optional adjustment)
  // arma::mat vhat;
  if (adjust) {
    return ((g.t() * (g)) / n) % ((c.t() * c) / (c.t() * c - 1));
  } else {
    return (g.t() * (g)) / n;
  }
}

arma::vec lambda2theta_ibd(const arma::vec& lambda,
                           const arma::vec& theta,
                           const arma::mat& g,
                           const arma::mat& c,
                           const double gamma) {
  // arma::vec arg = 1 + g * lambda;
  // arma::vec dplog_vec = dplog(arg, 1 / g.n_rows);
  arma::vec dplog_vec = dplog(1 + g * lambda, 1.0 / g.n_rows);
  // gradient
  arma::vec gradient = -arma::sum(arma::diagmat(dplog_vec) * c, 0).t() % lambda;
  // update theta by GD with lambda fixed
  // arma::vec theta_hat = theta - gamma * gradient;

  return theta - gamma * gradient;
}

arma::vec approx_lambda_ibd(const arma::mat& x,
                            const arma::mat& c,
                            const arma::vec& theta0,
                            const arma::vec& theta1,
                            const arma::vec& lambda0)
{
  arma::mat g0 = g_ibd(theta0, x, c);
  arma::vec arg = 1 + g0 * lambda0;

  // LHS
  arma::mat LHS = g0.t() * (g0.each_col() / arma::pow(arg, 2));

  // RHS
  arma::mat I_RHS = arma::diagmat(arma::sum((c.each_col() / arg), 0));
  arma::mat J = c.each_row() % arma::trans(lambda0);
  J.each_col() /= arma::pow(arg, 2);
  arma::mat J_RHS = g0.t() * J;
  arma::mat RHS = -I_RHS + J_RHS;

  // Jacobian matrix
  arma::mat jacobian = arma::solve(LHS, RHS, arma::solve_opts::fast);

  // linear approximation for lambda1
  return lambda0 + jacobian * (theta1 - theta0);
}

double cutoff_pairwise_PB_ibd(const arma::mat& x,
                              const arma::mat& c,
                              const int B,
                              const double level,
                              // const bool approx_lambda,
                              const bool adjust) {
  /// parameters ///
  const int p = x.n_cols;   // number of points(treatments)
  const std::vector<std::array<int, 2>> pairs = all_pairs(p);   // vector of pairs
  const int m = pairs.size();   // number of hypotheses
  const arma::mat V_hat = cov_ibd(x, c, adjust);    // covariance estimate

  /// A hat matrices ///
  arma::cube A_hat(p, p, m);
  for (int i = 0; i < m; ++i) {
    arma::rowvec R = arma::zeros(1, p);
    R(pairs[i][0] - 1) = 1;
    R(pairs[i][1] - 1) = -1;
    A_hat.slice(i) = (R.t() * R) / arma::as_scalar(R * V_hat * R.t());
  }

  // U hat matrices
  const arma::mat U_hat = arma::mvnrnd(arma::zeros(p), V_hat, B);

  // B bootstrap replicates(B x m matrix)
  arma::mat bootstrap_sample(B, m);
  for (int i = 0; i < m; ++i) {
    bootstrap_sample.col(i) =
      arma::diagvec(U_hat.t() * A_hat.slice(i) * U_hat);
  }

  return
    arma::as_scalar(arma::quantile(arma::max(bootstrap_sample, 1),
                                   arma::vec{1 - level}));
}

arma::mat centering_ibd(arma::mat x)
{
  // centering with nonzero elements
  x.each_col([](arma::vec& v) {
    v.elem(arma::find(v)) -= arma::mean(v.elem(arma::find(v)));
  });

  return x;
}

minEL test_ibd_EL(const arma::mat& x,
                  const arma::mat& c,
                  const arma::mat& L,
                  const arma::vec& rhs,
                  const bool approx_lambda,
                  const int maxit,
                  const double abstol) {
  /// initialization ///
  const int n = x.n_rows;
  const int r = arma::rank(L);
  if (r != L.n_rows) {
    Rcpp::stop("Hypothesis matrix L must have full rank.");
  }
  if (L.n_rows != rhs.n_elem) {
    Rcpp::stop("Dimensions of L and rhs do not match.");
  }
  // initial parameter value set as group means
  arma::vec theta = n * arma::trans(arma::mean(x, 0) / arma::sum(c, 0));
  // constraint imposed on the initial value by projection
  theta = linear_projection(theta, L, rhs);
  // estimating function
  arma::mat g = g_ibd(theta, x, c);
  // evaluation
  EL eval = getEL(g);
  arma::vec lambda = eval.lambda;
  // for current function value(-logLR)
  double f0 = arma::sum(plog(1 + g * lambda, 1.0 / n));
  // for updated function value
  double f1 = f0;

  /// minimization(projected gradient descent) ///
  double gamma = std::pow(arma::mean(arma::sum(c, 0)), -1);    // step size
  bool convergence = false;
  int iterations = 0;
  // proposed value for theta
  arma::vec theta_tmp;
  arma::vec lambda_tmp;
  arma::mat g_tmp;
  while (convergence == false) {
    if (f0 - f1 < abstol && iterations > 0) {
      convergence = true;
    } else {
      // update parameter by GD with lambda fixed
      theta_tmp = lambda2theta_ibd(lambda, theta, g, c, gamma);
      // projection
      theta_tmp = linear_projection(theta_tmp, L, rhs);
      // update g
      g_tmp = g_ibd(theta_tmp, x, c);
      if (approx_lambda && iterations > 1) {
        // update lambda
        lambda_tmp = approx_lambda_ibd(x, c, theta, theta_tmp, lambda);
      } else {
        // update lambda
        eval = getEL(g_tmp);
        lambda_tmp = eval.lambda;
        if (!eval.convergence && iterations > 9) {
          theta = theta_tmp;
          lambda = lambda_tmp;
          Rcpp::warning("Convex hull constraint not satisfied during optimization. Optimization halted.");
          break;
        }
      }
      // update function value
      f0 = f1;
      f1 = arma::sum(plog(1 + g_tmp * lambda_tmp, 1.0 / n));
      // step halving to ensure that the updated function value be
      // strinctly less than the current function value
      while (f0 <= f1) {
        // reduce step size
        gamma /= 2;
        // propose new theta
        theta_tmp = lambda2theta_ibd(lambda, theta, g, c, gamma);
        theta_tmp = linear_projection(theta_tmp, L, rhs);
        // propose new lambda
        g_tmp = g_ibd(theta_tmp, x, c);
        if (approx_lambda && iterations > 1) {
          lambda_tmp = approx_lambda_ibd(x, c, theta, theta_tmp, lambda);
        } else {
          eval = getEL(g_tmp);
          lambda_tmp = eval.lambda;
        }

        if (gamma < abstol) {
          // Rcpp::warning("Convex hull constraint not satisfied during step halving.");
          minEL result;
          result.theta = theta_tmp;
          result.lambda = lambda_tmp;
          result.nlogLR = f1;
          result.iterations = iterations;
          result.convergence = convergence;
          return result;
        }
        // propose new function value
        f1 = arma::sum(plog(1 + g_tmp * lambda_tmp, 1.0 / n));
      }
      // update parameters
      theta = theta_tmp;
      lambda = lambda_tmp;
      g = g_tmp;
      if (iterations == maxit) {
        break;
      }
      ++iterations;
    }
  }

  minEL result;
  result.theta = theta;
  result.lambda = lambda;
  result.nlogLR = f1;
  result.iterations = iterations;
  result.convergence = convergence;
  return result;
}

std::array<double, 2> pair_confidence_interval_ibd(const arma::mat& x,
                                                   const arma::mat& c,
                                                   const arma::mat& L,
                                                   const bool approx_lambda,
                                                   const double init,
                                                   const double threshold) {
  // upper endpoint
  double upper_lb = init;
  double upper_size = 0.5;
  double upper_ub = init + upper_size;
  double upper_eval =
    test_ibd_EL(x, c, L, arma::vec{upper_ub}, approx_lambda).nlogLR;
  // upper bound for upper endpoint
  while (2 * upper_eval <= threshold) {
    upper_size *= 2;
    upper_ub = init + upper_size;
    upper_eval =
      test_ibd_EL(x, c, L, arma::vec{upper_ub}, approx_lambda).nlogLR;
  }
  // approximate upper bound by numerical search
  while (upper_ub - upper_lb > 1e-04) {
    upper_eval =
      test_ibd_EL(x, c, L, arma::vec{(upper_lb + upper_ub) / 2}, approx_lambda).nlogLR;
    if (2 * upper_eval > threshold) {
      upper_ub = (upper_lb + upper_ub) / 2;
    } else {
      upper_lb = (upper_lb + upper_ub) / 2;
    }
  }

  // lower endpoint
  double lower_ub = init;
  double lower_size = 0.5;
  double lower_lb = init - lower_size;
  double lower_eval =
    test_ibd_EL(x, c, L, arma::vec{lower_lb}, approx_lambda).nlogLR;
  // lower bound for lower endpoint
  while (2 * lower_eval <= threshold) {
    lower_size *= 2;
    lower_lb = init - lower_size;
    lower_eval =
      test_ibd_EL(x, c, L, arma::vec{lower_lb}, approx_lambda).nlogLR;
  }
  // approximate lower bound by numerical search
  while (lower_ub - lower_lb > 1e-04) {
    lower_eval =
      test_ibd_EL(x, c, L, arma::vec{(lower_lb + lower_ub) / 2}, approx_lambda).nlogLR;
    if (2 * lower_eval > threshold) {
      lower_lb = (lower_lb + lower_ub) / 2;
    } else {
      lower_ub = (lower_lb + lower_ub) / 2;
    }
  }

  return std::array<double, 2>{lower_ub, upper_lb};
}

double cutoff_pairwise_NPB_ibd(const arma::mat& x,
                               const int B,
                               const double level,
                               const bool approx_lambda,
                               const int maxit,
                               const double abstol)
{
  // centered matrix
  arma::mat x_centered = centering_ibd(x);

  const int p = x.n_cols;
  const std::vector<std::array<int, 2>> pairs = all_pairs(p);   // vector of pairs
  const int m = pairs.size();   // number of hypotheses


  // B bootstrap test statistics(m x B matrix)
  arma::mat bootstrap_statistics(m, B);
  for (int b = 0; b < B; ++b) {
    arma::mat sample_b = bootstrap_sample(x_centered);
    arma::mat incidence_mat_b = arma::conv_to<arma::mat>::from(sample_b != 0);
    for (int j = 0; j < m; ++j) {
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[j][0] - 1) = 1;
      L(pairs[j][1] - 1) = -1;
      const minEL pairwise_result =
        test_ibd_EL(sample_b, incidence_mat_b, L, arma::zeros(1), approx_lambda, maxit, abstol);
      // We do not check the convex hull constraint when NPB is used.
      // if (!pairwise_result.convergence) {
      //   Rcpp::warning("Test for pair (%i,%i) failed in bootstrap sample %i. \n",
      //                 pairs[j][0], pairs[j][1], b + 1);
      // }
      bootstrap_statistics(j, b) = 2 * pairwise_result.nlogLR;
    }
    if (b % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  return
    arma::as_scalar(arma::quantile(arma::max(bootstrap_statistics, 0),
                                   arma::vec{1 - level}));
}
