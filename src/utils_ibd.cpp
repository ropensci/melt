#include "utils_ibd.h"

Eigen::MatrixXd g_ibd(const Eigen::Ref<const Eigen::VectorXd>& theta,
                      const Eigen::Ref<const Eigen::MatrixXd>& x,
                      const Eigen::Ref<const Eigen::MatrixXd>& c) {
  return x - (c.array().rowwise() * theta.array().transpose()).matrix();
}

Eigen::MatrixXd cov_ibd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                        const Eigen::Ref<const Eigen::MatrixXd>& c) {
  // // estimator(global minimizer)
  // const Eigen::VectorXd theta_hat =
  //   x.array().colwise().sum() / c.array().colwise().sum();
  // estimating function
  Eigen::MatrixXd g =
    g_ibd(x.array().colwise().sum() / c.array().colwise().sum(), x, c);
  // covariance estimate
  return (g.transpose() * g) / x.rows();
}

Eigen::VectorXd lambda2theta_ibd(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::VectorXd>& theta,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const double gamma) {
  // Eigen::VectorXd dplog_vec =
  //   PSEUDO_LOG::dp(Eigen::VectorXd::Ones(g.rows()) + g * lambda);
  // // gradient
  // Eigen::VectorXd gradient =
  //   -(dplog_vec.asDiagonal() * c).array().colwise().sum().transpose() * lambda.array();
  // // update theta by GD with lambda fixed
  // return theta - gamma * gradient;

  Eigen::VectorXd ngradient =
    (PSEUDO_LOG::dp(Eigen::VectorXd::Ones(g.rows()) + g * lambda).matrix().asDiagonal() * c)
  .array().colwise().sum().transpose() * lambda.array();
  return theta + gamma * ngradient;
}

void lambda2theta_void(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    Eigen::Ref<Eigen::VectorXd> theta,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const double gamma) {
  // Eigen::VectorXd ngradient =
  //   (PSEUDO_LOG::dp(
  //       Eigen::VectorXd::Ones(g.rows()) + g * lambda).matrix().asDiagonal() * c)
  //                                 .array().colwise().sum().transpose() * lambda.array();
  theta +=
    gamma * ((PSEUDO_LOG::dp(
    Eigen::VectorXd::Ones(
      g.rows()) + g * lambda).matrix().asDiagonal() * c)
               .array().colwise().sum().transpose() *
      lambda.array()).matrix();
}

Eigen::VectorXd approx_lambda_ibd(
    const Eigen::Ref<const Eigen::MatrixXd>& g0,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const Eigen::Ref<const Eigen::VectorXd>& theta0,
    const Eigen::Ref<const Eigen::VectorXd>& theta1,
    const Eigen::Ref<const Eigen::VectorXd>& lambda0) {
  Eigen::ArrayXd&& arg = Eigen::VectorXd::Ones(g0.rows()) + g0 * lambda0;
  Eigen::ArrayXd&& denominator = Eigen::pow(arg, 2);

  // LHS
  Eigen::MatrixXd&& LHS =
    g0.transpose() * (g0.array().colwise() / denominator).matrix();

  // RHS
  const Eigen::MatrixXd I_RHS =
    ((c.array().colwise() / arg).colwise().sum()).matrix().asDiagonal();
  const Eigen::MatrixXd J_RHS =
    (g0.array().colwise() / denominator).matrix().transpose() *
    (c.array().rowwise() * lambda0.array().transpose()).matrix();
  Eigen::MatrixXd&& RHS = -I_RHS + J_RHS;

  // Jacobian matrix
  Eigen::MatrixXd&& jacobian = LHS.ldlt().solve(RHS);

  // linear approximation for lambda1
  return lambda0 + jacobian * (theta1 - theta0);
}

std::array<double, 2> pair_confidence_interval_ibd(
    const Eigen::Ref<const Eigen::VectorXd>& theta0,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const double init,
    const double threshold) {
  // upper endpoint
  double upper_lb = init;
  double upper_size = 1;
  double upper_ub = init + upper_size;
  // upper bound for upper endpoint
  while (2 * test_nlogLR(theta0, x, c,
                         lhs, Eigen::Matrix<double, 1, 1>(upper_ub)) <= threshold) {
    upper_lb = upper_ub;
    upper_ub += upper_size;
  }
  // approximate upper bound by numerical search
  while (upper_ub - upper_lb > 1e-04) {
    if (2 * test_nlogLR(theta0, x, c, lhs,
                        Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2)) > threshold) {
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
  while (2 * test_nlogLR(theta0, x, c,
                         lhs, Eigen::Matrix<double, 1, 1>(lower_lb)) <= threshold) {
    lower_ub = lower_lb;
    lower_lb -= lower_size;
  }
  // approximate lower bound by numerical search
  while (lower_ub - lower_lb > 1e-04) {
    if (2 * test_nlogLR(theta0, x, c, lhs,
                        Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2)) > threshold) {
      lower_lb = (lower_lb + lower_ub) / 2;
    } else {
      lower_ub = (lower_lb + lower_ub) / 2;
    }
  }
  return std::array<double, 2>{lower_ub, upper_lb};
}

Eigen::MatrixXd centering_ibd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                              const Eigen::Ref<const Eigen::MatrixXd>& c) {
  return x -
    (c.array().rowwise() *
    (x.array().colwise().sum() / c.array().colwise().sum())).matrix();
}

Eigen::MatrixXd rmvn(const Eigen::MatrixXd& x, const int n) {
  // generate standard multivariate gaussian random vectors(n by p matrix)
  Eigen::MatrixXd I(n, x.cols());
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < n; ++i) {
      I(i, j) = R::rnorm(0, 1.0);
    }
  }
  // get the square root matrix of the covariance matrix
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(x);
  // return the target normal random vectors(n by p matrix)
  return I * es.operatorSqrt();
}

Eigen::ArrayXd bootstrap_statistics_pairwise_AMC(
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const int k,
    const std::vector<std::array<int, 2>>& pairs,
    const int B,
    const double level) {
  // covariance estimate
  const Eigen::MatrixXd V_hat = cov_ibd(x, c);
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

double cutoff_pairwise_NB_approx(const Eigen::Ref<const Eigen::MatrixXd>& x,
                                 const Eigen::Ref<const Eigen::MatrixXd>& c,
                                 const int k,
                                 const std::vector<std::array<int, 2>>& pairs,
                                 const int B,
                                 const double level,
                                 const int ncores,
                                 const int maxit,
                                 const double abstol) {
  const int n = x.rows();
  const int p = x.cols();
  const int m = pairs.size();   // number of hypotheses

  // centered matrix
  // Eigen::MatrixXd&& x_centered = centering_ibd(x, c);
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

  // B bootstrap results(we only need maximum statistics)
  Eigen::VectorXd bootstrap_statistics(B);
  #pragma omp parallel for num_threads(ncores) default(none) shared(B, maxit, abstol, pairs, x_centered, c, p, m, bootstrap_index, bootstrap_statistics) schedule(auto)
  for (int b = 0; b < B; ++b) {
    Eigen::ArrayXd statistics_b(m);
    for (int j = 0; j < m; ++j) {
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
      lhs(pairs[j][0]) = 1;
      lhs(pairs[j][1]) = -1;
      statistics_b(j) =
        2 * test_ibd_EL_approx(
            bootstrap_sample(x_centered, bootstrap_index.col(b)),
            bootstrap_sample(c, bootstrap_index.col(b)),
            lhs, Eigen::Matrix<double, 1, 1>(0),
            maxit, abstol).nlogLR;
      }
    // need to generalize later for k-FWER control
    bootstrap_statistics(b) = statistics_b.maxCoeff();
    }

  // quantile function needed!
  Rcpp::Function quantile("quantile");
  return
    Rcpp::as<double>(quantile(bootstrap_statistics,
                              Rcpp::Named("probs") = 1 - level));
}

Eigen::ArrayXd bootstrap_statistics_pairwise_NB(
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const int k,
    const std::vector<std::array<int, 2>>& pairs,
    const int B,
    const double level,
    const int ncores,
    const int maxit,
    const double abstol) {
  const int n = x.rows();
  const int p = x.cols();
  const int m = pairs.size();   // number of hypotheses

  // centered matrix
  // Eigen::MatrixXd&& x_centered = centering_ibd(x, c);
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
  Eigen::ArrayXd k_bootstrap_statistics(B);
  #pragma omp parallel for num_threads(ncores) default(none) shared(x_centered, c, k, pairs, B, maxit, abstol, p, m, bootstrap_index, k_bootstrap_statistics) schedule(auto)
  for (int b = 0; b < B; ++b) {
    std::vector<double> bootstrap_statistics(m);
    for (int j = 0; j < m; ++j) {
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
      lhs(pairs[j][0]) = 1;
      lhs(pairs[j][1]) = -1;
      bootstrap_statistics[j] =
        2 * test_nlogLR(bootstrap_sample(x_centered, bootstrap_index.col(b)),
                        bootstrap_sample(c, bootstrap_index.col(b)),
                        lhs, Eigen::Matrix<double, 1, 1>(0),
                        maxit, abstol);
    }
    // kth largest element
    std::sort(bootstrap_statistics.begin(), bootstrap_statistics.end());
    k_bootstrap_statistics[b] = bootstrap_statistics[m - k];
  }
  return k_bootstrap_statistics;
}

minEL test_ibd_EL(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                  const Eigen::Ref<const Eigen::MatrixXd>& x,
                  const Eigen::Ref<const Eigen::MatrixXd>& c,
                  const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                  const Eigen::Ref<const Eigen::VectorXd>& rhs,
                  const int maxit,
                  const double abstol) {
  /// initialization ///
  // Constraint imposed on the initial value by projection.
  // The initial value is given as treatment means.
  Eigen::VectorXd theta =
    linear_projection(theta0, lhs, rhs);
  // estimating function
  Eigen::MatrixXd g = g_ibd(theta, x, c);
  // evaluation
  Eigen::VectorXd lambda = EL2(g).lambda;
  // for current function value(-logLR)
  double f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * lambda);

  /// minimization(projected gradient descent) ///
  double gamma = 1.0 / (c.colwise().sum().mean());    // step size
  bool convergence = false;
  int iterations = 0;
  // proposed value for theta
  while (!convergence && iterations != maxit) {
    // update parameter by GD with lambda fixed -> projection
    Eigen::VectorXd theta_tmp = theta;
    lambda2theta_void(lambda, theta_tmp, g, c, gamma);
    linear_projection_void(theta_tmp, lhs, rhs);
    // update g
    Eigen::MatrixXd g_tmp = g_ibd(theta_tmp, x, c);
    // update lambda
    EL2 eval(g_tmp);
    Eigen::VectorXd lambda_tmp = eval.lambda;
    if (!eval.convergence && iterations > 9) {
      lambda = std::move(lambda_tmp);
      Rcpp::warning("Convex hull constraint not satisfied during optimization. Optimization halted.");
      break;
    }

    // update function value
    double f0 = f1;
    f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);

    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < f1) {
      // reduce step size
      gamma /= 2;
      // propose new theta
      theta_tmp = theta;
      lambda2theta_void(lambda, theta_tmp, g, c, gamma);
      linear_projection_void(theta_tmp, lhs, rhs);
      // propose new lambda
      g_tmp = g_ibd(theta_tmp, x, c);
      lambda_tmp = EL2(g_tmp).lambda;
      if (gamma < abstol) {
        lambda = std::move(lambda_tmp);
        Rcpp::warning("Convex hull constraint not satisfied during step halving.");
        break;
      }
      // propose new function value
      f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    }

    // update parameters
    theta = std::move(theta_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);
    ++iterations;

    // convergence check
    if (f0 - f1 < abstol) {
      convergence = true;
    }
  }

  return {theta, lambda, f1, iterations, convergence};
}

double test_nlogLR(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                   const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const int maxit,
                   const double abstol) {
  /// initialization ///
  // Constraint imposed on the initial value by projection.
  // The initial value is given as treatment means.
  Eigen::VectorXd theta =
    linear_projection(theta0, lhs, rhs);
  // estimating function
  Eigen::MatrixXd g = g_ibd(theta, x, c);
  // evaluation
  Eigen::VectorXd lambda = EL2(g).lambda;
  // for current function value(-logLR)
  double f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * lambda);

  /// minimization(projected gradient descent) ///
  double gamma = 1.0 / (c.colwise().sum().mean());    // step size
  bool convergence = false;
  int iterations = 0;
  // proposed value for theta
  while (!convergence && iterations != maxit) {
    // update parameter by GD with lambda fixed -> projection
    Eigen::VectorXd theta_tmp = theta;
    lambda2theta_void(lambda, theta_tmp, g, c, gamma);
    linear_projection_void(theta_tmp, lhs, rhs);
    // update g
    Eigen::MatrixXd g_tmp = g_ibd(theta_tmp, x, c);
    // update lambda
    EL2 eval(g_tmp);
    Eigen::VectorXd lambda_tmp = eval.lambda;
    if (!eval.convergence && iterations > 9) {
      lambda = std::move(lambda_tmp);
      Rcpp::warning("Convex hull constraint not satisfied during optimization. Optimization halted.");
      break;
    }

    // update function value
    double f0 = f1;
    f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);

    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < f1) {
      // reduce step size
      gamma /= 2;
      // propose new theta
      theta_tmp = theta;
      lambda2theta_void(lambda, theta_tmp, g, c, gamma);
      linear_projection_void(theta_tmp, lhs, rhs);
      // propose new lambda
      g_tmp = g_ibd(theta_tmp, x, c);
      lambda_tmp = EL2(g_tmp).lambda;
      if (gamma < abstol) {
        lambda = std::move(lambda_tmp);
        Rcpp::warning("Convex hull constraint not satisfied during step halving.");
        break;
      }
      // propose new function value
      f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    }

    // update parameters
    theta = std::move(theta_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);
    ++iterations;

    // convergence check
    if (f0 - f1 < abstol) {
      convergence = true;
    }
  }

  return f1;
}

double test_nlogLR(const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const int maxit,
                   const double abstol) {
  /// initialization ///
  // Constraint imposed on the initial value by projection.
  // The initial value is given as treatment means.
  Eigen::VectorXd theta =
    linear_projection(x.array().colwise().sum() / c.array().colwise().sum(),
                      lhs, rhs);
  // estimating function
  Eigen::MatrixXd g = g_ibd(theta, x, c);
  // evaluation
  Eigen::VectorXd lambda = EL2(g).lambda;
  // for current function value(-logLR)
  double f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * lambda);

  /// minimization(projected gradient descent) ///
  double gamma = 1.0 / (c.colwise().sum().mean());    // step size
  bool convergence = false;
  int iterations = 0;
  // proposed value for theta
  while (!convergence && iterations != maxit) {
    // update parameter by GD with lambda fixed -> projection
    Eigen::VectorXd theta_tmp = theta;
    lambda2theta_void(lambda, theta_tmp, g, c, gamma);
    linear_projection_void(theta_tmp, lhs, rhs);
    // update g
    Eigen::MatrixXd g_tmp = g_ibd(theta_tmp, x, c);
    // update lambda
    EL2 eval(g_tmp);
    Eigen::VectorXd lambda_tmp = eval.lambda;
    if (!eval.convergence && iterations > 9) {
      lambda = std::move(lambda_tmp);
      Rcpp::warning("Convex hull constraint not satisfied during optimization. Optimization halted.");
      break;
    }
    // update function value
    double f0 = f1;
    f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < f1) {
      // reduce step size
      gamma /= 2;
      // propose new theta
      theta_tmp = theta;
      lambda2theta_void(lambda, theta_tmp, g, c, gamma);
      linear_projection_void(theta_tmp, lhs, rhs);
      // propose new lambda
      g_tmp = g_ibd(theta_tmp, x, c);
      lambda_tmp = EL2(g_tmp).lambda;
      if (gamma < abstol) {
        lambda = std::move(lambda_tmp);
        Rcpp::warning("Convex hull constraint not satisfied during step halving.");
        break;
      }
      // propose new function value
      f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    }

    // update parameters
    theta = std::move(theta_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);
    ++iterations;

    // convergence check
    if (f0 - f1 < abstol) {
      convergence = true;
    }
  }

  return f1;
}

//
minEL test_ibd_EL_approx(const Eigen::Ref<const Eigen::MatrixXd>& x,
                         const Eigen::Ref<const Eigen::MatrixXd>& c,
                         const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                         const Eigen::Ref<const Eigen::VectorXd>& rhs,
                         const int maxit,
                         const double abstol) {
  /// initialization ///
  // Constraint imposed on the initial value by projection.
  // The initial value is given as treatment means.
  Eigen::VectorXd theta =
    linear_projection(x.array().colwise().sum() / c.array().colwise().sum(),
                      lhs, rhs);

  // estimating function
  Eigen::MatrixXd g = g_ibd(theta, x, c);
  // evaluation
  EL eval = getEL(g);
  Eigen::VectorXd lambda = eval.lambda;
  // for current function value(-logLR)
  double f0 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * lambda);
  // for updated function value
  double f1 = f0;

  /// minimization(projected gradient descent) ///
  double gamma = 1.0 / (c.array().colwise().sum().mean());    // step size
  bool convergence = false;
  int iterations = 0;

  while (!convergence && iterations != maxit) {
    // update parameter by GD with lambda fixed -> projection
    Eigen::VectorXd theta_tmp =
      linear_projection(lambda2theta_ibd(lambda, theta, g, c, gamma), lhs, rhs);
    // update g
    Eigen::MatrixXd g_tmp = g_ibd(theta_tmp, x, c);

    Eigen::VectorXd lambda_tmp(theta.size());
      if (iterations > 1) {
      // update lambda
      lambda_tmp = approx_lambda_ibd(g, c, theta, theta_tmp, lambda);
    } else {
      // update lambda
      eval = getEL(g_tmp);
      lambda_tmp = eval.lambda;
      if (!eval.convergence && iterations > 9) {
        theta = std::move(theta_tmp);
        lambda = std::move(lambda_tmp);
        Rcpp::warning("Convex hull constraint not satisfied during optimization. Optimization halted.");
        break;
      }
    }

    // update function value
    f0 = f1;
    f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);

    // step halving to ensure that the updated function value be
    // strictly less than the current function value
    while (f0 < f1) {
      // reduce step size
      gamma /= 2;
      // propose new theta
      theta_tmp =
        linear_projection(lambda2theta_ibd(lambda, theta, g, c, gamma),
                          lhs, rhs);
      // propose new lambda
      g_tmp = g_ibd(theta_tmp, x, c);
      if (iterations > 1) {
        lambda_tmp = approx_lambda_ibd(g, c, theta, theta_tmp, lambda);
      } else {
        eval = getEL(g_tmp);
        lambda_tmp = eval.lambda;
      }
      if (gamma < abstol) {
        theta = std::move(theta_tmp);
        lambda = std::move(lambda_tmp);
        Rcpp::warning("Convex hull constraint not satisfied during step halving.");
        break;
      }
      // propose new function value
      f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    }

    // update parameters
    theta = std::move(theta_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);

    // convergence check
    if (f0 - f1 < abstol && iterations > 0) {
      convergence = true;
    } else {
      ++iterations;
    }
  }

  return {theta, lambda, f1, iterations, convergence};
}
