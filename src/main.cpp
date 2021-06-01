#include "utils.h"
#include "utils_ibd.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


//' Compute empirical likelihood for mean
//'
//' Compute empirical likelihood for mean
//'
//' @param theta a vector of parameters to be tested.
//' @param x a matrix or vector of data. Each row is an observation vector.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 50.
//' @export
// [[Rcpp::export]]
Rcpp::List el_mean(arma::rowvec theta, arma::mat x,
               int maxit = 100, double abstol = 1e-8) {
  if (arma::rank(x) != x.n_cols) {
    Rcpp::stop("Design matrix x must have full rank.");
  }
  // ncol must be same as p
  int n = x.n_rows;
  int p = theta.n_elem;

  // estimating equation
  arma::mat g = x;
  g.each_row() -= theta;

  // minimization
  arma::vec l; l.zeros(p);
  arma::vec lc;
  arma::vec arg = 1 + g * l;
  arma::vec y;
  arma::mat J;
  arma::mat Q;
  arma::mat R;
  double f0;
  double f1;
  int iterations = 0;
  bool convergence = false;
  while (convergence == false) {
    // function evaluation(initial)
    f0 = -arma::sum(plog(arg, 1 / n));
    // J matrix & y vector
    arma::vec v1 = arma::sqrt(-d2plog(arg, 1 / n));
    arma::vec v2 = dplog(arg, 1 / n);
    J = g.each_col() % v1;
    y = v2 / v1;
    // update lambda by NR method with least square & step halving
    arma::qr_econ(Q, R, J);
    lc = l + arma::solve(R, trans(Q) * y);
    double alpha = 1;
    while(-arma::sum(plog(1 + g * lc, 1 / n)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * solve(R, trans(Q) * y);
    }
    // update function value
    l = lc;
    arg = 1 + g * l;
    f1 = -arma::sum(plog(arg, 1 / n));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      arma::vec v1 = arma::sqrt(-d2plog(arg, 1 / n));
      arma::vec v2 = dplog(arg, 1 / n);
      J = g.each_col() % v1;
      y = v2 / v1;
      convergence = true;
    } else {
      ++iterations;
      if(iterations == maxit) {
        break;
      }
    }
  }

  // result
  return Rcpp::List::create(
    Rcpp::Named("nlogLR") = -f1,
    Rcpp::Named("lambda") = Rcpp::NumericVector(l.begin(), l.end()),
    Rcpp::Named("grad") = Rcpp::as<std::vector<double>>(Rcpp::wrap(-trans(J) * y)),
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("convergence") = convergence);
}

//' Empirical likelihood test for mean
//'
//' Compute empirical likelihood for mean
//'
//' @param theta a vector of parameters to be tested.
//' @param x a matrix or vector of data. Each row is an observation vector.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 50.
//' @export
// [[Rcpp::export]]
Rcpp::List el_mean2(const arma::vec& theta,
                    const arma::mat& x,
                    int maxit = 100,
                    double abstol = 1e-8) {
  if (arma::rank(x) != x.n_cols) {
    Rcpp::stop("Design matrix x must have full rank.");
  }
  // estimating function for mean parameters
  arma::mat g = g_mean(theta, x);

  // compute EL
  EL result = getEL(g, maxit, abstol);

  return Rcpp::List::create(
    Rcpp::Named("nlogLR") = result.nlogLR,
    Rcpp::Named("lambda") =
      Rcpp::NumericVector(result.lambda.begin(), result.lambda.end()),
    Rcpp::Named("gradient") =
      Rcpp::NumericVector(result.gradient.begin(), result.gradient.end()),
    Rcpp::Named("iterations") = result.iterations,
    Rcpp::Named("convergence") = result.convergence);
}

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
                    const double& abstol = 1e-8) {
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
                        const bool approx_lambda = false,
                        const int maxit = 1e4,
                        const double abstol = 1e-8) {
  if (level <= 0 || level >= 1)
  {
    Rcpp::stop("level must be between 0 and 1.");
  }
  if (method != "PB" && method != "NPB")
  {
    Rcpp::warning
    ("method '%s' is not supported. Using 'PB' as default.",
     method);
    method = "PB";
  }

  const int n = x.n_rows;
  const int p = x.n_cols;
  // all pairs
  std::vector<std::vector<int>> pairs = all_pairs(p);
  // number of hypotheses
  const int m = pairs.size();
  // global minimizer
  arma::vec theta_hat = n * arma::trans(arma::mean(x, 0) / arma::sum(c, 0));
  // estimate
  Rcpp::NumericVector estimate(m);
  // statistics(-2logLR)
  Rcpp::NumericVector statistic(m);
  // cutoff value
  double cutoff;
  if (method == "NPB")
  {
    cutoff = cutoff_pairwise_NPB_ibd(x, B, level, approx_lambda, maxit, abstol);
  }
  else
  {
    cutoff = cutoff_pairwise_PB_ibd(x, c, B, level, vcov_adj);
  }

  if (!interval)
  {
    for(int i = 0; i < m; ++i)
    {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      minEL pairwise_result =
        test_ibd_EL(x, c, L, arma::zeros(1), false, maxit, abstol);
      if (!pairwise_result.convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      double nlogLR = pairwise_result.nlogLR;
      statistic(i) = 2 * pairwise_result.nlogLR;
    }

    Rcpp::List model_info =
      Rcpp::List::create(Rcpp::Named("model.matrix") = x ,
                         Rcpp::Named("incidence.matrix") = c);
    Rcpp::List result =
      Rcpp::List::create(
        Rcpp::Named("estimate") = estimate,
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("level") = level,
        Rcpp::Named("cutoff") = cutoff,
        Rcpp::Named("method") = method,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") = model_info);
    result.attr("class") = "pairwise.ibd";

    return result;
  }
  else
  {
    Rcpp::List CI(m);
    for(int i = 0; i < m; ++i)
    {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      minEL pairwise_result =
        test_ibd_EL(x, c, L, arma::zeros(1), false, maxit, abstol);
      if (!pairwise_result.convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      double nlogLR = pairwise_result.nlogLR;
      statistic(i) = 2 * pairwise_result.nlogLR;

      // confidence interval(optional)
      CI(i) =
        pair_confidence_interval_ibd(x, c, L, approx_lambda, estimate(i), cutoff);
    }
    Rcpp::List model_info =
      Rcpp::List::create(Rcpp::Named("model.matrix") = x ,
                         Rcpp::Named("incidence.matrix") = c);
    Rcpp::List result =
      Rcpp::List::create(
        Rcpp::Named("estimate") = estimate,
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("CI") = CI,
        Rcpp::Named("level") = level,
        Rcpp::Named("cutoff") = cutoff,
        Rcpp::Named("method") = method,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") = model_info);
    result.attr("class") = "pairwise.ibd";

    return result;
  }
}





//
Rcpp::List pairwise_PB_ibd(const arma::mat& x,
                           const arma::mat& c,
                           const bool& interval = false,
                           const int B = 1e5,
                           const double& level = 0.05,
                           const bool vcov_adj = false,
                           const bool approx_lambda = false,
                           const int maxit = 1e3,
                           const double& abstol = 1e-8) {
  const int n = x.n_rows;
  const int p = x.n_cols;
  // all pairs
  std::vector<std::vector<int>> pairs = all_pairs(p);
  // number of hypotheses
  const int m = pairs.size();
  // global minimizer
  arma::vec theta_hat = n * arma::trans(arma::mean(x, 0) / arma::sum(c, 0));
  // estimate
  Rcpp::NumericVector estimate(m);
  // statistics(-2logLR)
  Rcpp::NumericVector statistic(m);

  if (interval) {
    if (level <= 0 || level >= 1) {
      Rcpp::stop("level must be between 0 and 1.");
    }
    double cutoff =
      cutoff_pairwise_PB_ibd(x, c, B, level, vcov_adj);
    Rcpp::List CI(m);
    for(int i = 0; i < m; ++i) {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      minEL pairwise_result =
        test_ibd_EL(x, c, L, arma::zeros(1), maxit, abstol);
      bool convergence = pairwise_result.convergence;
      if (!convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      double nlogLR = pairwise_result.nlogLR;
      statistic(i) = 2 * nlogLR;

      // CI(optional)
      if (interval) {
        CI(i) =
          pair_confidence_interval_ibd(x, c, L, approx_lambda, estimate(i), cutoff);
      }
    }
    Rcpp::List model_info =
      Rcpp::List::create(Rcpp::Named("model.matrix") = x ,
                         Rcpp::Named("incidence.matrix") = c);
    Rcpp::List result =
      Rcpp::List::create(
        Rcpp::Named("estimate") = estimate,
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("CI") = CI,
        Rcpp::Named("level") = level,
        Rcpp::Named("cutoff") = cutoff,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") = model_info);
    result.attr("class") = "pairwise.ibd";
    return result;

  } else {
    double cutoff = cutoff_pairwise_PB_ibd(x, c, B, level, vcov_adj);
    for(int i = 0; i < m; ++i) {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      minEL pairwise_result =
        test_ibd_EL(x, c, L, arma::zeros(1), maxit, abstol);
      bool convergence = pairwise_result.convergence;
      if (!convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      double nlogLR = pairwise_result.nlogLR;
      statistic(i) = 2 * nlogLR;
    }
    Rcpp::List model_info =
      Rcpp::List::create(Rcpp::Named("model.matrix") = x ,
                         Rcpp::Named("incidence.matrix") = c);
    Rcpp::List result =
      Rcpp::List::create(
        Rcpp::Named("estimate") = estimate,
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("level") = level,
        Rcpp::Named("cutoff") = cutoff,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") = model_info);
    result.attr("class") = "pairwise.ibd";
    return result;
  }
}

Rcpp::List minP_pairwise_ibd(const arma::mat& x,
                        const arma::mat& c,
                        const bool interval = false,
                        const int B = 1e4,
                        const double level = 0.05,
                        const bool approx_lambda = false,
                        const int maxit = 1e4,
                        const double abstol = 1e-8) {
  if (level <= 0 || level >= 1)
  {
    Rcpp::stop("level must be between 0 and 1.");
  }

  const int n = x.n_rows;
  const int p = x.n_cols;
  // all pairs
  std::vector<std::vector<int>> pairs = all_pairs(p);
  // number of hypotheses
  const int m = pairs.size();
  // global minimizer
  arma::vec theta_hat = n * arma::trans(arma::mean(x, 0) / arma::sum(c, 0));
  // estimate
  Rcpp::NumericVector estimate(m);
  // statistics(-2logLR)
  Rcpp::NumericVector statistic(m);
  // p-values(F-calibrated)
  Rcpp::NumericVector p_value(m);
  // common quantile value
  double quantile =
    quantile_pairwise_NPB_ibd(x, B, level, maxit, abstol);


  if (!interval)
  {
    for(int i = 0; i < m; ++i)
    {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      minEL pairwise_result =
        test_ibd_EL(x, c, L, arma::zeros(1), maxit, abstol);
      if (!pairwise_result.convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      statistic(i) = 2 * pairwise_result.nlogLR;
      p_value(i) = pairwise_result.p_value;
    }

    Rcpp::List model_info =
      Rcpp::List::create(Rcpp::Named("model.matrix") = x ,
                         Rcpp::Named("incidence.matrix") = c);
    Rcpp::List result =
      Rcpp::List::create(
        Rcpp::Named("estimate") = estimate,
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("p.val") = p_value,
        Rcpp::Named("quantile") = quantile,
        Rcpp::Named("level") = level,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") = model_info);
    result.attr("class") = "pairwise.ibd";

    return result;
  }
  else
  {
    Rcpp::List CI(m);
    for(int i = 0; i < m; ++i)
    {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      minEL pairwise_result =
        test_ibd_EL(x, c, L, arma::zeros(1), maxit, abstol);
      if (!pairwise_result.convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      double nlogLR = pairwise_result.nlogLR;
      statistic(i) = 2 * pairwise_result.nlogLR;

      // confidence interval(optional)
      // CI(i) = pair_confidence_interval_ibd(x, c, L, estimate(i), cutoff);
    }
    Rcpp::List model_info =
      Rcpp::List::create(Rcpp::Named("model.matrix") = x ,
                         Rcpp::Named("incidence.matrix") = c);
    Rcpp::List result =
      Rcpp::List::create(
        Rcpp::Named("estimate") = estimate,
        Rcpp::Named("statistic") = statistic,
        Rcpp::Named("CI") = CI,
        Rcpp::Named("level") = level,
        // Rcpp::Named("cutoff") = cutoff,
        // Rcpp::Named("method") = method,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") = model_info);
    result.attr("class") = "pairwise.ibd";

    return result;
  }
}

arma::mat tt(const arma::mat& x,
                                 const int B,
                                 const double level,
                                 const int maxit,
                                 const double abstol)
{
  // centered matrix(null transformation)
  arma::mat x_centered = centering_ibd(x);
  const int n = x.n_rows;
  const int p = x.n_cols;
  const std::vector<std::vector<int>> pairs = all_pairs(p);   // vector of pairs
  const int m = pairs.size();   // number of hypotheses

  // B bootstrap test statistics(B x m matrix)
  arma::mat bootstrap_statistics(B, m);
  for (int b = 0; b < B; ++b) {
    arma::mat sample_b = bootstrap_sample(x_centered);
    arma::mat incidence_mat_b = arma::conv_to<arma::mat>::from(sample_b != 0);
    for (int j = 0; j < m; ++j)
    {
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[j][0] - 1) = 1;
      L(pairs[j][1] - 1) = -1;
      const minEL& pairwise_result =
        test_ibd_EL(sample_b, incidence_mat_b, L, arma::zeros(1), maxit, abstol);
      if (!pairwise_result.convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed in %i bootstrap sample. \n",
                      pairs[j][0], pairs[j][1], b);
      }
      bootstrap_statistics(b, j) = 2 * pairwise_result.nlogLR;
    }
    if (b % 100 == 0)
    {
      Rcpp::checkUserInterrupt();
    }
  }

  // bootstrap unadjusted p-values based on rank of statistics
  arma::mat bootstrap_pvalue_unadjusted(B, m);
  for (int j = 0; j < m; ++j)
  {
    bootstrap_pvalue_unadjusted.col(j) =
      (arma::conv_to<arma::vec>::from(
          arma::sort_index(
            arma::sort_index(bootstrap_statistics.col(j), "descend"),
            "ascend") + 1)) / B;
    // note that we add 1
  }
  // // F-calibration(df1 = 1, df2 = n - 1, lower = false, log = false)
  //
  // bootstrap_pvalue_unadjusted(b, j) =
  //   R::pf(2.0 * pairwise_result.nlogLR, 1, n - 1, false, false);
  // return
  // arma::as_scalar(bootstrap_pvalue_unadjusted(0,0));
  return bootstrap_pvalue_unadjusted;
  // return
  // arma::as_scalar(arma::quantile(arma::min(bootstrap_pvalue_unadjusted, 1),
  // arma::vec{level}));
}
