#include "utils.h"

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
      iterations++;
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
                    const int& maxit = 100,
                    const double& abstol = 1e-8) {
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
                    const int& maxit = 1000,
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
      iterations++;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("theta") = Rcpp::NumericVector(theta.begin(), theta.end()),
    Rcpp::Named("lambda") = Rcpp::NumericVector(lambda.begin(), lambda.end()),
    Rcpp::Named("nlogLR") = f1,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("convergence") = convergence);
}

std::vector<double> pair_confidence_interval(const arma::mat& x,
                                             const arma::mat& c,
                                             const arma::mat& L,
                                             const double& init,
                                             const double& threshold) {
  // upper endpoint
  double upper_lb = init;
  double upper_size = 1;
  double upper_ub = init + upper_size;
  double upper_eval = test_ibd(x, c, L, arma::vec{upper_ub})["nlogLR"];
  // upper bound for upper endpoint
  while (2 * upper_eval <= threshold) {
    upper_size *= 2;
    upper_ub = init + upper_size;
    upper_eval = test_ibd(x, c, L, arma::vec{upper_ub})["nlogLR"];
  }
  // approximate upper bound by numerical search
  while (upper_ub - upper_lb >= 1e-04) {
    upper_eval =
      test_ibd(x, c, L, arma::vec {(upper_lb + upper_ub) / 2})["nlogLR"];
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
  double lower_eval = test_ibd(x, c, L, arma::vec{lower_lb})["nlogLR"];
  // lower bound for lower endpoint
  while (2 * lower_eval <= threshold) {
    lower_size *= 2;
    lower_lb = init - lower_size;
    lower_eval = test_ibd(x, c, L, arma::vec{lower_lb})["nlogLR"];
  }
  // approximate lower bound by numerical search
  while (lower_ub - lower_lb >= 1e-04) {
    lower_eval =
      test_ibd(x, c, L, arma::vec {(lower_lb + lower_ub) / 2})["nlogLR"];
    if (2 * lower_eval <= threshold) {
      lower_ub = (lower_lb + lower_ub) / 2;
    } else {
      lower_lb = (lower_lb + lower_ub) / 2;
    }
  }

  std::vector<double> confidence_inverval{lower_ub, upper_lb};
  return confidence_inverval;
}

//' Pairwise comparison for BIBD
//'
//' Pairwise comparison for BIBD
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param interval whether to compute interval.
//' @param B number of bootstrap replicate.
//' @param level confidence level.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List pairwise_ibd(const arma::mat& x,
                        const arma::mat& c,
                        const bool& interval = false,
                        const int& B = 1e5,
                        const double& level = 0.95,
                        const int& maxit = 1e3,
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
    if (level < 0 || level > 1) {
      Rcpp::stop("level must be between 0 and 1.");
    }
    double threshold = threshold_pairwise_ibd(x, c, B, level);
    Rcpp::List CI(m);
    for(int i = 0; i < m; i++) {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      Rcpp::List pairwise_result =
        test_ibd(x, c, L, arma::zeros(1), maxit, abstol);
      bool convergence = pairwise_result["convergence"];
      if (!convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      double nlogLR = pairwise_result["nlogLR"];
      statistic(i) = 2 * nlogLR;

      // CI(optional)
      if (interval) {
        CI(i) = pair_confidence_interval(x, c, L, estimate(i), threshold);
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
        Rcpp::Named("threshold") = threshold,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") = model_info);
    result.attr("class") = "pairwise.ibd";
    return result;

  } else {
    double threshold = threshold_pairwise_ibd(x, c, B, level);
    for(int i = 0; i < m; i++) {
      // estimates
      estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);

      // statistics
      arma::rowvec L = arma::zeros(1, p);
      L(pairs[i][0] - 1) = 1;
      L(pairs[i][1] - 1) = -1;
      Rcpp::List pairwise_result =
        test_ibd(x, c, L, arma::zeros(1), maxit, abstol);
      bool convergence = pairwise_result["convergence"];
      if (!convergence) {
        Rcpp::warning("Test for pair (%i,%i) failed. \n",
                      pairs[i][0], pairs[i][1]);
      }
      double nlogLR = pairwise_result["nlogLR"];
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
        Rcpp::Named("threshold") = threshold,
        Rcpp::Named("num.bootstrap") = B,
        Rcpp::Named("model.info") = model_info);
    result.attr("class") = "pairwise.ibd";
    return result;
  }
}
