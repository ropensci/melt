#include "utils_ibd.h"

//' Hypothesis test for incomplete block design
//'
//' Hypothesis test for incomplete block design
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param L a linear hypothesis matrix.
//' @param rhs right-hand-side vector for hypothesis, with as many entries as rows in the hypothesis matrix.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//' @export
// [[Rcpp::export]]
Rcpp::List test_ibd(const arma::mat& x,
                    const arma::mat& c,
                    const arma::mat& L,
                    const arma::vec& rhs,
                    const int maxit = 1000,
                    const double abstol = 1e-8) {
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
          return Rcpp::List::create(
            Rcpp::Named("theta") = Rcpp::NumericVector(theta.begin(), theta.end()),
            Rcpp::Named("lambda") = Rcpp::NumericVector(lambda.begin(), lambda.end()),
            Rcpp::Named("nlogLR") = f1,
            Rcpp::Named("iterations") = iterations,
            Rcpp::Named("convergence") = convergence);
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

  return Rcpp::List::create(
    Rcpp::Named("theta") = Rcpp::NumericVector(theta.begin(), theta.end()),
    Rcpp::Named("lambda") = Rcpp::NumericVector(lambda.begin(), lambda.end()),
    Rcpp::Named("nlogLR") = f1,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("convergence") = convergence);
}

//' Pairwise comparison for Incomplete Block Design
//'
//' Pairwise comparison for Incomplete Block Design
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param interval whether to compute interval. Defaults to FALSE.
//' @param B number of bootstrap replicates.
//' @param level level.
//' @param method the method to be used; either 'PB' or 'NPB' is supported. Defaults to 'PB'.
//' @param vcov_adj whether to adjust for the covariance estimate. Defaults to FALSE.
//' @param approx_lambda whether to use the approximation for lambda. Defaults to FALSE.
//' @param ncores number of cores(threads) to use. Defaults to 1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List pairwise_ibd(const arma::mat& x,
                        const arma::mat& c,
                        const bool interval = false,
                        const int B = 1e4,
                        const double level = 0.05,
                        std::string method = "PB",
                        const bool vcov_adj = false,
                        const bool cheat = false,
                        const bool approx_lambda = false,
                        const int ncores = 1,
                        const int maxit = 1e4,
                        const double abstol = 1e-8) {
  if (level <= 0 || level >= 1) {
    Rcpp::stop("level must be between 0 and 1.");
  }
  if (method != "PB" && method != "NPB") {
    Rcpp::warning
    ("method '%s' is not supported. Using 'PB' as default.",
     method);
    method = "PB";
  }

  const int n = x.n_rows;
  const int p = x.n_cols;
  // all pairs
  std::vector<std::array<int, 2>> pairs = all_pairs(p);
  // number of hypotheses
  const int m = pairs.size();
  // global minimizer
  const arma::vec theta_hat =
    n * arma::trans(arma::mean(x, 0) / arma::sum(c, 0));
  // estimate
  Rcpp::NumericVector estimate(m);
  // statistics(-2logLR)
  Rcpp::NumericVector statistic(m);
  // cutoff value
  const double cutoff =
    method == "NPB" ?
    cutoff_pairwise_NPB_ibd(x, B, level, approx_lambda, ncores, maxit, abstol) :
    cutoff_pairwise_PB_ibd(x, c, B, level, vcov_adj, cheat);

  if (!interval) {
    for (int i = 0; i < m; ++i) {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      const minEL pairwise_result =
        test_ibd_EL(x, c, L, arma::zeros(1), false, maxit, abstol);
      if (!pairwise_result.convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      statistic(i) = 2 * pairwise_result.nlogLR;
    }

    Rcpp::List result =
      Rcpp::List::create(
        Rcpp::Named("estimate") = estimate,
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("level") = level,
        Rcpp::Named("cutoff") = cutoff,
        Rcpp::Named("method") = method,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") =
          Rcpp::List::create(Rcpp::Named("model.matrix") = x,
                             Rcpp::Named("incidence.matrix") = c));
    result.attr("class") = "pairwise.ibd";

    return result;
  } else {
    Rcpp::List CI(m);
    for (int i = 0; i < m; ++i) {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      const minEL pairwise_result =
        test_ibd_EL(x, c, L, arma::zeros(1), false, maxit, abstol);
      if (!pairwise_result.convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      statistic(i) = 2 * pairwise_result.nlogLR;

      // confidence interval(optional)
      CI(i) =
        pair_confidence_interval_ibd(
          x, c, L, approx_lambda, estimate(i), cutoff);
    }

    Rcpp::List result =
      Rcpp::List::create(
        Rcpp::Named("estimate") = estimate,
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("CI") = CI,
        Rcpp::Named("level") = level,
        Rcpp::Named("cutoff") = cutoff,
        Rcpp::Named("method") = method,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") =
          Rcpp::List::create(Rcpp::Named("model.matrix") = x ,
                             Rcpp::Named("incidence.matrix") = c));
    result.attr("class") = "pairwise.ibd";

    return result;
  }
}
