#include "utils_ibd.h"

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat g_ibd(const arma::vec& theta,
                const arma::mat& x,
                const arma::mat& c) {
  return x - c.each_row() % theta.t();
}


arma::mat cov_ibd(const arma::mat& x,
                  const arma::mat& c,
                  const bool adjust) {
  // number of blocks
  int n = x.n_rows;
  // estimator(global minimizer)
  arma::vec theta_hat = n * arma::trans(arma::mean(x, 0) / arma::sum(c, 0));
  // estimating function
  arma::mat g = x - c.each_row() % theta_hat.t();
  // covariance estimate(optional adjustment)
  arma::mat vhat;
  if (adjust) {
    vhat = ((g.t() * (g)) / n) % ((c.t() * c) / (c.t() * c - 1));
  } else {
    vhat = (g.t() * (g)) / n;
  }

  return vhat;
}

arma::vec lambda2theta_ibd(const arma::vec& lambda,
                           const arma::vec& theta,
                           const arma::mat& g,
                           const arma::mat& c,
                           const double& gamma) {
  arma::vec arg = 1 + g * lambda;
  arma::vec dplog_vec = dplog(arg, 1 / g.n_rows);
  // gradient
  arma::vec gradient = -arma::sum(arma::diagmat(dplog_vec) * c, 0).t() % lambda;
  // update theta by GD with lambda fixed
  arma::vec theta_hat = theta - gamma * gradient;

  return theta_hat;
}

double cutoff_pairwise_PB_ibd(const arma::mat& x,
                              const arma::mat& c,
                              const int B,
                              const double level,
                              const bool adjust) {
  /// parameters ///
  const int p = x.n_cols;   // number of points(treatments)
  const std::vector<std::vector<int>> pairs = all_pairs(p);   // vector of pairs
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
                  const int maxit,
                  const double abstol) {
  /// initialization ///
  const int n = x.n_rows;
  if (arma::rank(L) != L.n_rows) {
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
  // If the convex hull constraint is not satisfied at the initial value, end.
  arma::vec arg = 1 + g * lambda;
  // for current function value(-logLR)
  double f0 = arma::sum(plog(arg, 1 / n));
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
      // update lambda
      eval = getEL(g_tmp);
      lambda_tmp = eval.lambda;
      if (!eval.convergence && iterations > 9) {
        theta = theta_tmp;
        lambda = lambda_tmp;
        Rcpp::warning("Convex hull constraint not satisfied during optimization. Optimization halted.");
        break;
      }
      // update function value
      f0 = f1;
      arg = 1 + g_tmp * lambda_tmp;
      f1 = arma::sum(plog(arg, 1 / n));
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
        eval = getEL(g_tmp);
        lambda_tmp = eval.lambda;
        if (gamma < abstol) {
          theta = theta_tmp;
          lambda = lambda_tmp;
          Rcpp::warning("Convex hull constraint not satisfied during step halving.");

          minEL result;
          result.theta = theta;
          result.lambda = lambda;
          result.nlogLR = f1;
          result.iterations = iterations;
          result.convergence = convergence;
          return result;
        }
        // propose new function value
        arg = 1 + g_tmp * lambda_tmp;
        f1 = arma::sum(plog(arg, 1 / n));
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

std::vector<double> pair_confidence_interval_ibd(const arma::mat& x,
                                                 const arma::mat& c,
                                                 const arma::mat& L,
                                                 const double& init,
                                                 const double& threshold) {
  // upper endpoint
  double upper_lb = init;
  double upper_size = 1;
  double upper_ub = init + upper_size;
  double upper_eval = test_ibd_EL(x, c, L, arma::vec{upper_ub}).nlogLR;
  // upper bound for upper endpoint
  while (2 * upper_eval <= threshold) {
    upper_size *= 2;
    upper_ub = init + upper_size;
    upper_eval = test_ibd_EL(x, c, L, arma::vec{upper_ub}).nlogLR;
  }
  // approximate upper bound by numerical search
  while (upper_ub - upper_lb >= 1e-04) {
    upper_eval =
      test_ibd_EL(x, c, L, arma::vec {(upper_lb + upper_ub) / 2}).nlogLR;
    if (2 * upper_eval <= threshold) {
      upper_lb = (upper_lb + upper_ub) / 2;
    } else {
      upper_ub = (upper_lb + upper_ub) / 2;
    }
  }

  // lower endpoint
  double lower_ub = init;
  double lower_size = 1;
  double lower_lb = init - lower_size;
  double lower_eval = test_ibd_EL(x, c, L, arma::vec{lower_lb}).nlogLR;
  // lower bound for lower endpoint
  while (2 * lower_eval <= threshold) {
    lower_size *= 2;
    lower_lb = init - lower_size;
    lower_eval = test_ibd_EL(x, c, L, arma::vec{lower_lb}).nlogLR;
  }
  // approximate lower bound by numerical search
  while (lower_ub - lower_lb >= 1e-04) {
    lower_eval =
      test_ibd_EL(x, c, L, arma::vec {(lower_lb + lower_ub) / 2}).nlogLR;
    if (2 * lower_eval <= threshold) {
      lower_ub = (lower_lb + lower_ub) / 2;
    } else {
      lower_lb = (lower_lb + lower_ub) / 2;
    }
  }

  std::vector<double> confidence_inverval{lower_ub, upper_lb};
  return confidence_inverval;
}

double cutoff_pairwise_NPB_ibd(const arma::mat& x,
                               const int B,
                               const double level,
                               const int maxit,
                               const double abstol)
{
  // centered matrix
  arma::mat x_centered = centering_ibd(x);

  const int p = x.n_cols;
  const std::vector<std::vector<int>> pairs = all_pairs(p);   // vector of pairs
  const int m = pairs.size();   // number of hypotheses


  // B bootstrap test statistics(B x m matrix)
  arma::mat bootstrap_statistics(B, m);
  for (int b = 0; b < B; ++b) {
    arma::mat sample_b = bootstrap_sample(x_centered);
    arma::mat incidence_mat_b = arma::conv_to<arma::mat>::from(sample_b != 0);
    for(int j = 0; j < m; ++j)
    {
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[j][0] - 1) = 1;
      L(pairs[j][1] - 1) = -1;
      minEL pairwise_result =
        test_ibd_EL(sample_b, incidence_mat_b, L, arma::zeros(1), maxit, abstol);
      if (!pairwise_result.convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed in %i bootstrap sample. \n",
                      pairs[j][0], pairs[j][1], b);
      }
      bootstrap_statistics(b, j) = 2 * pairwise_result.nlogLR;
    }
  }

  return
    arma::as_scalar(arma::quantile(arma::max(bootstrap_statistics, 1),
                                   arma::vec{1 - level}));
}
