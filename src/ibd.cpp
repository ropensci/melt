#include "utils_ibd.h"

//' Hypothesis test for incomplete block design
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
Rcpp::List test_ibd(const Eigen::MatrixXd& x,
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
    test_ibd_EL(x.array().colwise().sum() / c.array().colwise().sum(),
                x, c, lhs, rhs, maxit, abstol);

  return Rcpp::List::create(
    Rcpp::Named("theta") = result.theta,
    Rcpp::Named("lambda") = result.lambda,
    Rcpp::Named("nlogLR") = result.nlogLR,
    Rcpp::Named("iterations") = result.iterations,
    Rcpp::Named("convergence") = result.convergence);
}

//' Pairwise comparison for Incomplete Block Design
//'
//' Pairwise comparison for Incomplete Block Design
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param interval whether to compute interval. Defaults to FALSE.
//' @param nbootstrap number of bootstrap replicates.
//' @param level level.
//' @param method the method to be used; either 'PB' or 'NB' is supported. Defaults to 'PB'.
//' @param correction whether to use blocked bootstrap. Defaults to FALSE.
//' @param approx whether to use the approximation for lambda. Defaults to FALSE.
//' @param nthread number of cores(threads) to use. Defaults to 1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List pairwise_ibd(const Eigen::MatrixXd& x,
                        const Eigen::MatrixXd& c,
                        const bool interval = false,
                        const int nbootstrap = 1e4,
                        const double level = 0.05,
                        std::string method = "PB",
                        const bool correction = false,
                        const bool approx = false,
                        const int nthread = 1,
                        const int maxit = 1e4,
                        const double abstol = 1e-8) {
  if (level <= 0 || level >= 1) {
    Rcpp::stop("level must be between 0 and 1.");
  }
  // all pairs
  std::vector<std::array<int, 2>> pairs = all_pairs(x.cols());
  // cutoff value
  double cutoff;
  if (method == "PB") {
    cutoff = cutoff_pairwise_PB(x, c, pairs, nbootstrap, level, correction);
  } else if (method != "NB") {
    Rcpp::warning
    ("method '%s' is not supported. Using 'PB' as default.",
     method);
    method = "PB";
    cutoff = cutoff_pairwise_PB(x, c, pairs, nbootstrap, level, correction);
  } else if (approx) {
    cutoff = cutoff_pairwise_NB_approx(x, c, nbootstrap, level, nthread, maxit, abstol);
  } else {
    cutoff = cutoff_pairwise_NB(x, c, nbootstrap, level, nthread, maxit, abstol);
  }
  // global minimizer
  const Eigen::VectorXd theta_hat =
    x.array().colwise().sum() / c.array().colwise().sum();
  // number of hypotheses
  const int m = pairs.size();


  // estimates
  Rcpp::NumericVector estimate(m);

  // statistics(-2logLR)
  Rcpp::NumericVector statistic(m);
  for (int i = 0; i < m; ++i) {
    // estimates
    estimate(i) = theta_hat(pairs[i][0] - 1) - theta_hat(pairs[i][1] - 1);
    // statistics
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
    lhs(pairs[i][0] - 1) = 1;
    lhs(pairs[i][1] - 1) = -1;
    minEL pairwise_result =
      test_ibd_EL(theta_hat, x, c, lhs, Eigen::Matrix<double, 1, 1>(0),
                  maxit, abstol);
    if (!pairwise_result.convergence) {
      Rcpp::warning("test for pair (%i,%i) failed. \n",
                    pairs[i][0], pairs[i][1]);
    }
    statistic(i) = 2 * pairwise_result.nlogLR;
  }

  // // cutoff value
  // double cutoff;
  // if (method == "PB") {
  //   cutoff = cutoff_pairwise_PB(x, c, pairs, nbootstrap, level, correction);
  // } else if (approx) {
  //   cutoff = cutoff_pairwise_NB_approx(x, c, nbootstrap, level, nthread, maxit, abstol);
  // } else {
  //   cutoff = cutoff_pairwise_NB(x, c, nbootstrap, level, nthread, maxit, abstol);
  // }

  // result
  Rcpp::List result;
  result["estimate"] = estimate;
  result["statistic"] = statistic;
  // confidence interval(optional)
  if (interval) {
    Rcpp::List CI(m);
    for (int i = 0; i < m; ++i) {
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
      lhs(pairs[i][0] - 1) = 1;
      lhs(pairs[i][1] - 1) = -1;
      CI(i) =
        pair_confidence_interval_ibd(theta_hat, x, c, lhs, estimate(i), cutoff);
    }
    result["CI"] = CI;
  }
  result["level"] = level;
  result["cutoff"] = cutoff;
  result["method"] = method;
  // result["num.bootstrap"] = nbootstrap;
  result.attr("class") = "elmulttest";
  return result;
}

//
//
//
//' Pairwise comparison for Incomplete Block Design
//'
//' Pairwise comparison for Incomplete Block Design
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param interval whether to compute interval. Defaults to FALSE.
//' @param B number of bootstrap replicates.
//' @param level level.
//' @param method the method to be used; either 'PB' or 'NB' is supported. Defaults to 'PB'.
//' @param correction whether to use blocked bootstrap. Defaults to FALSE.
//' @param approx whether to use the approximation for lambda. Defaults to FALSE.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List tt(const Eigen::MatrixXd& x,
                        const Eigen::MatrixXd& c,
                        const int maxit = 1e4,
                        const double abstol = 1e-8) {
  // all pairs
  std::vector<std::array<int, 2>> pairs = all_pairs(x.cols());

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
    lhs(pairs[i][0] - 1) = 1;
    lhs(pairs[i][1] - 1) = -1;
    minEL pairwise_result =
      test_ibd_EL(theta_hat, x, c, lhs, Eigen::Matrix<double, 1, 1>(0),
                  maxit, abstol);
    if (!pairwise_result.convergence) {
      Rcpp::warning("Test for pair (%i,%i) failed. \n",
                    pairs[i][0], pairs[i][1]);
    }
    statistic(i) = 2 * pairwise_result.nlogLR;
  }

  // result
  Rcpp::List result;
  result["statistic"] = statistic;
  result.attr("class") = "pairwise.ibd";
  return result;
}
