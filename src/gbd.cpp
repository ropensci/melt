#include "utils_gbd.h"

//' Empirical Likelihood Hypothesis Testing
//'
//' Hypothesis test for incomplete block design
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param lhs a linear hypothesis matrix.
//' @param rhs right-hand-side vector for hypothesis, with as many entries as rows in the hypothesis matrix.
//' @param approx whether to use the approximation for lambda. Defaults to FALSE.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//' @export
// [[Rcpp::export]]
Rcpp::List test_gbd(const Eigen::MatrixXd& x,
                    const Eigen::MatrixXd& c,
                    const Eigen::MatrixXd& lhs,
                    const Eigen::VectorXd& rhs,
                    const bool approx = false,
                    const int maxit = 1000,
                    const double abstol = 1e-8) {
  /// initialization ///
  // if (arma::rank(L) != lhs.rows()) {
  //   Rcpp::stop("Hypothesis matrix lhs must have full rank.");
  // }
  if (lhs.rows() != rhs.rows()) {
    Rcpp::stop("Dimensions of L and rhs do not match.");
  }
  minEL result =
    test_gbd_EL(x.array().colwise().sum() / c.array().colwise().sum(),
                x, c, lhs, rhs, maxit, abstol);

  return Rcpp::List::create(
    Rcpp::Named("theta") = result.theta,
    Rcpp::Named("lambda") = result.lambda,
    Rcpp::Named("nlogLR") = result.nlogLR,
    Rcpp::Named("iterations") = result.iterations,
    Rcpp::Named("convergence") = result.convergence);
}

//' Pairwise Comparisons for General Block Design
//'
//' Either all pairwise comparisons or comparisons with control is available.
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param control control treatment. Defaults to 0.
//' @param k integer k for k-FWER. Defaults to 1.
//' @param interval whether to compute interval. Defaults to TRUE.
//' @param B number of bootstrap replicates.
//' @param level level.
//' @param method the method to be used; either 'AMC' or 'NB' is supported. Defaults to 'AMC'.
//' @param approx whether to use the approximation for lambda. Defaults to FALSE.
//' @param nthread number of cores(threads) to use. Defaults to 1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List el_pairwise(const Eigen::MatrixXd& x,
                       const Eigen::MatrixXd& c,
                       const int control = 0,
                       const int k = 1,
                       const double level = 0.05,
                       const bool interval = true,
                       const std::string method = "AMC",
                       const int B = 1e4,
                       const bool approx = false,
                       const int nthread = 1,
                       const int maxit = 1e4,
                       const double abstol = 1e-8) {
  if (level <= 0 || level >= 1) {
    Rcpp::stop("level must be between 0 and 1.");
  }
  // pairs
  std::vector<std::array<int, 2>> pairs = comparison_pairs(x.cols(), control);
  // number of hypotheses
  const int m = pairs.size();

  // bootstrap statitics
  Eigen::ArrayXd bootstrap_statistics_pairwise(B);
  if (method == "AMC") {
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_AMC(x, c, k, pairs, B, level);
  } else if (approx) {
    // NOT READY
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_NB(x, c, k, pairs, B,
                                       level, nthread, maxit, abstol);
  } else {
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_NB(x, c, k, pairs, B,
                                       level, nthread, maxit, abstol);
  }

  // cutoff
  Rcpp::Function quantile("quantile");
  const double cutoff =
    Rcpp::as<double>(quantile(bootstrap_statistics_pairwise,
                              Rcpp::Named("probs") = 1 - level));

  // estimates
  std::vector<double> estimate(m);

  // statistics
  std::vector<double> statistic(m);

  // convergences
  std::vector<double> convergence(m);

  // adjusted p-values
  std::vector<double> adj_pvalues(m);

  // global minimizer
  const Eigen::VectorXd theta_hat =
    x.array().colwise().sum() / c.array().colwise().sum();
  for (int i = 0; i < m; ++i) {
    estimate[i] = theta_hat(pairs[i][0]) - theta_hat(pairs[i][1]);

    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
    lhs(pairs[i][0]) = 1;
    lhs(pairs[i][1]) = -1;
    minEL pairwise_result =
      test_gbd_EL(theta_hat, x, c, lhs, Eigen::Matrix<double, 1, 1>(0),
                  maxit, abstol);
    if (!pairwise_result.convergence) {
      Rcpp::warning("test for pair (%i,%i) failed. \n",
                    pairs[i][0] + 1, pairs[i][1] + 1);
    }
    statistic[i] = 2 * pairwise_result.nlogLR;

    convergence[i] = pairwise_result.convergence;

    adj_pvalues[i] =
      static_cast<double>(
        (bootstrap_statistics_pairwise >= statistic[i]).count()) / B;
  }

  // result
  Rcpp::List result;
  result["estimate"] = estimate;
  result["statistic"] = statistic;
  result["convergence"] = convergence;
  // confidence interval(optional)
  if (interval) {
    std::vector<double> lower(m);
    std::vector<double> upper(m);
    for (int i = 0; i < m; ++i) {
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
      lhs(pairs[i][0]) = 1;
      lhs(pairs[i][1]) = -1;
      std::array<double, 2> ci =
        pair_confidence_interval_gbd(theta_hat, x, c, lhs, estimate[i], cutoff);
      lower[i] = ci[0];
      upper[i] = ci[1];
    }
    result["lower"] = lower;
    result["upper"] = upper;
  }
  result["p.adj"] = adj_pvalues;
  result["k"] = k;
  result["level"] = level;
  result["method"] = approx ? "NB(approx)" : method;
  result["cutoff"] = cutoff;
  result.attr("class") = Rcpp::CharacterVector({"pairwise", "elmulttest"});
  return result;
}

// [[Rcpp::export]]
Rcpp::List tt(const Eigen::MatrixXd& x,
              const Eigen::MatrixXd& c,
              const int control = 0,
              const int k = 1,
              const int maxit = 1e4,
              const double abstol = 1e-8) {
  // all pairs
  std::vector<std::array<int, 2>> pairs = comparison_pairs(x.cols(), control);

  // global minimizer
  const Eigen::VectorXd theta_hat =
    x.array().colwise().sum() / c.array().colwise().sum();
  // number of hypotheses
  const int m = pairs.size();

  // statistics(-2logLR)
  Rcpp::NumericVector statistic(m);
  for (int i = 0; i < m; ++i) {
    // statistics
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
    lhs(pairs[i][0]) = 1;
    lhs(pairs[i][1]) = -1;
    minEL pairwise_result =
      test_gbd_EL(theta_hat, x, c, lhs, Eigen::Matrix<double, 1, 1>(0),
                  maxit, abstol);
    if (!pairwise_result.convergence) {
      Rcpp::warning("Test for pair (%i,%i) failed. \n",
                    pairs[i][0] + 1, pairs[i][1] + 1);
    }
    statistic(i) = 2 * pairwise_result.nlogLR;
  }

  // result
  Rcpp::List result;
  result["statistic"] = statistic;
  result.attr("class") = "pairwise.gbd";
  return result;
}
