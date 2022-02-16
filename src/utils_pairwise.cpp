#include "utils_pairwise.h"

std::vector<std::array<int, 2>> comparison_pairs(
    const int p, const int control) {
  // initialize a vector of vectors
  std::vector<std::array<int, 2>> pairs;
  if (control == 0){
    // the size of vector is p choose 2
    pairs.reserve(p * (p - 1) / 2);
    // fill in each elements(pairs)
    for (int i = 0; i < p - 1; ++i) {
      for (int j = i + 1; j < p; ++j) {
        pairs.emplace_back(std::array<int, 2>{i, j});

      }
    }
  } else {
    // the size of vector is p - 1
    pairs.reserve(p - 1);
    // fill in each elements(pairs)
    for (int i = 0; i < p; ++i) {
      if (i == control - 1) {
        continue;
      }
      pairs.emplace_back(std::array<int, 2>{i, control - 1});
    }
  }
  return pairs;
}

std::array<double, 2> pair_confidence_interval_gbd(
    const Eigen::Ref<const Eigen::VectorXd>& theta0,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const double threshold,
    const double init,
    const double cutoff) {
  // upper endpoint
  double upper_lb = init;
  double upper_size = 1;
  double upper_ub = init + upper_size;
  // upper bound for upper endpoint
  while (2 * test_nlogLR(
      theta0, x, c, lhs,
      Eigen::Matrix<double, 1, 1>(upper_ub), 1000, 1e-04, threshold) <= cutoff) {
    upper_lb = upper_ub;
    upper_ub += upper_size;
  }
  // approximate upper bound by numerical search
  while (upper_ub - upper_lb > 1e-04) {
    if (2 * test_nlogLR(
        theta0, x, c, lhs,
        Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2), 1000, 1e-04,
        threshold) > cutoff) {
      upper_ub = (upper_lb + upper_ub) / 2;
    } else {
      upper_lb = (upper_lb + upper_ub) / 2;
    }
  }

  // // lower endpoint
  // double lower_ub;
  // double lower_size = upper_ub - init;
  // double lower_lb = init - lower_size;

  // lower endpoint
  double lower_ub = init;
  double lower_size = 1;
  double lower_lb = init - lower_size;


  // lower bound for lower endpoint
  while (2 * test_nlogLR(
      theta0, x, c, lhs,
      Eigen::Matrix<double, 1, 1>(lower_lb), 1000, 1e-04, threshold) <= cutoff) {
    lower_ub = lower_lb;
    lower_lb -= lower_size;
  }
  // approximate lower bound by numerical search
  while (lower_ub - lower_lb > 1e-04) {
    if (2 * test_nlogLR(
        theta0, x, c, lhs,
        Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2), 1000, 1e-04,
        threshold) > cutoff) {
      lower_lb = (lower_lb + lower_ub) / 2;
    } else {
      lower_ub = (lower_lb + lower_ub) / 2;
    }
  }
  return std::array<double, 2>{lower_ub, upper_lb};
}

Eigen::ArrayXd bootstrap_statistics_pairwise_AMC(
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const int k,
    const std::vector<std::array<int, 2>>& pairs,
    const int B,
    const double level) {
  // covariance estimate
  const Eigen::MatrixXd V_hat = cov_gbd(x, c);
  // U hat matrices
  const Eigen::MatrixXd U_hat = rmvn(V_hat, B);
  // number of hypotheses
  const int m = pairs.size();

  // B bootstrap statistics(B x m matrix)
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                Eigen::RowMajor> bootstrap_statistics(B, m);
  for (int j = 0; j < m; ++j) {
    Eigen::RowVectorXd R = Eigen::RowVectorXd::Zero(1, x.cols());
    R(pairs[j][0]) = 1;
    R(pairs[j][1]) = -1;
    Eigen::MatrixXd A_hat = (R.transpose() * R) / (R * V_hat * R.transpose());
    bootstrap_statistics.col(j) =
      (U_hat * A_hat * U_hat.transpose()).diagonal();
  }
  // rowwise sort
  for (int b = 0; b < B; ++b) {
    std::sort(bootstrap_statistics.row(b).data(),
              bootstrap_statistics.row(b).data() + m);
  }
  return bootstrap_statistics.col(m - k);
}

Eigen::ArrayXd bootstrap_statistics_pairwise_NB(
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const int k,
    const std::vector<std::array<int, 2>>& pairs,
    const int B,
    const double level,
    const int nthread,
    const bool progress,
    const double threshold,
    const int maxit,
    const double abstol) {
  const int n = x.rows();
  const int p = x.cols();
  const int m = pairs.size();   // number of hypotheses

  // centered matrix
  // Eigen::MatrixXd&& x_centered = centering_gbd(x, c);
  const Eigen::MatrixXd x_centered =
    x - (c.array().rowwise() *
    (x.array().colwise().sum() / c.array().colwise().sum())).matrix();

  // index vector for boostrap(length n * B)
  // generate index to sample(Rcpp) -> transform to std::vector ->
  // reshape to ArrayXXi(Eigen)
  const Eigen::ArrayXXi bootstrap_index =
    Eigen::Map<Eigen::ArrayXXi, Eigen::Unaligned>(
        (Rcpp::as<std::vector<int>>(
            Rcpp::sample(Rcpp::IntegerVector(Rcpp::seq(0, n - 1)), n * B, true)))
    .data(), n, B);

  // B bootstrap results
  Eigen::ArrayXd k_bootstrap_statistics = Eigen::ArrayXd::Constant(B, NA_REAL);
  NB_ProgressBar pb;
  Progress pg(B, progress, pb);
  bool stop = false;
  #pragma omp parallel for num_threads(nthread) schedule(dynamic)
  for (int b = 0; b < B; ++b) {
    if (!pg.is_aborted()) { // the only way to exit an OpenMP loop
      // check for user abort for every 250 iterations
      // if ((b + 1) % 250 == 0){
      if (Progress::check_abort()) {
          stop = true;
      }
      // }

      std::vector<double> bootstrap_statistics(m);
      for (int j = 0; j < m; ++j) {
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
      lhs(pairs[j][0]) = 1;
      lhs(pairs[j][1]) = -1;
      bootstrap_statistics[j] =
        2 * test_nlogLR(bootstrap_sample(x_centered, bootstrap_index.col(b)),
                        bootstrap_sample(c, bootstrap_index.col(b)),
                        lhs, Eigen::Matrix<double, 1, 1>(0),
                        threshold, maxit, abstol);
        }
    // kth largest element
    std::sort(bootstrap_statistics.begin(), bootstrap_statistics.end());
    k_bootstrap_statistics[b] = bootstrap_statistics[m - k];
    // update the progress
    pg.increment();
    }
  }
  if (stop) {
    Rcpp::NumericVector v_interrupted = Rcpp::wrap(k_bootstrap_statistics);
    v_interrupted = v_interrupted[!Rcpp::is_na(v_interrupted)];
    k_bootstrap_statistics = Rcpp::as<Eigen::ArrayXd>(v_interrupted);
    if (progress) {
      REprintf("\ninterrupted!");
    }
  }
  return k_bootstrap_statistics;
}
