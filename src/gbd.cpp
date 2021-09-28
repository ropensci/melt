#include "utils_gbd.h"

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
                       const bool progress = true,
                       const double threshold = 50,
                       const int maxit = 1e4,
                       const double abstol = 1e-8) {
  if (level <= 0 || level >= 1) {
    Rcpp::stop("level must be between 0 and 1.");
  }
  // pairs
  std::vector<std::array<int, 2>> pairs = comparison_pairs(x.cols(), control);
  // number of hypotheses
  const int m = pairs.size();
  // estimates
  std::vector<double> estimate(m);
  // statistics
  std::vector<double> statistic(m);
  // convergences
  std::vector<bool> convergence(m);

  // global minimizer
  const Eigen::VectorXd theta_hat =
    x.array().colwise().sum() / c.array().colwise().sum();
  if (progress) {
    REprintf("computing statistics...");
  }
  // statistics
  for (int i = 0; i < m; ++i) {
    Rcpp::checkUserInterrupt();
    estimate[i] = theta_hat(pairs[i][0]) - theta_hat(pairs[i][1]);
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
    lhs(pairs[i][0]) = 1;
    lhs(pairs[i][1]) = -1;
    minEL pairwise_result =
      test_gbd_EL(theta_hat, x, c, lhs, Eigen::Matrix<double, 1, 1>(0),
                  threshold, maxit, abstol);
    statistic[i] = 2 * pairwise_result.nlogLR;
    convergence[i] = pairwise_result.convergence;
  }

  // bootstrap statitics
  if (progress) {
    REprintf("\ncomputing cutoff...");
  }
  Eigen::ArrayXd bootstrap_statistics_pairwise(B);
  if (method == "AMC") {
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_AMC(x, c, k, pairs, B, level);
  } else {
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_NB(x, c, k, pairs, B,
                                       level, nthread, progress,
                                       threshold, maxit, abstol);
  }
  // adjusted p-values
  std::vector<double> adj_pvalues(m);
  for (int i = 0; i < m; ++i) {
    adj_pvalues[i] =
      static_cast<double>(
        (bootstrap_statistics_pairwise >= statistic[i]).count()) / B;
  }
  // cutoff
  bool rr = std::any_of(convergence.begin(), convergence.end(),
                        [](bool v) {return !v;});
  if (rr) {
    // Rcpp::Rcout << "sdf\n";
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_AMC(x, c, k, pairs, B, level);
  }
  Rcpp::Function quantile("quantile");
  const double cutoff =
    Rcpp::as<double>(quantile(bootstrap_statistics_pairwise,
                              Rcpp::Named("probs") = 1 - level));

  // result
  Rcpp::List result;
  result["estimate"] = estimate;
  result["statistic"] = statistic;
  result["convergence"] = convergence;
  // confidence interval(optional)
  if (progress) {
    REprintf("\ncomputing confidence intervals...\n");
  }
  if (interval) {
    std::vector<double> lower(m);
    std::vector<double> upper(m);
    for (int i = 0; i < m; ++i) {
      Rcpp::checkUserInterrupt();
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
      lhs(pairs[i][0]) = 1;
      lhs(pairs[i][1]) = -1;
      std::array<double, 2> ci =
        pair_confidence_interval_gbd(theta_hat, x, c, lhs,
                                     threshold, estimate[i], cutoff);
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
  result["B"] = bootstrap_statistics_pairwise.size();
  result["cutoff"] = cutoff;
  result.attr("class") = Rcpp::CharacterVector({"pairwise", "melt"});
  return result;
}
