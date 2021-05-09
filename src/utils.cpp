#include <RcppArmadillo.h>
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::NumericVector plog(Rcpp::NumericVector x, double threshold) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = std::log(x[i]);
    } else {
      out[i] = std::log(threshold) - 1.5 + 2 * std::pow(threshold, -1) * x[i] -
        std::pow(x[i] / threshold, 2) / 2;
    }
  }
  return out;
}

Rcpp::NumericVector dplog(Rcpp::NumericVector x, double threshold) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = std::pow(x[i], -1);
    } else {
      out[i] = 2 * std::pow(threshold, -1) - x[i] * std::pow(threshold, -2);
    }
  }
  return out;
}

Rcpp::NumericVector d2plog(Rcpp::NumericVector x, double threshold) {
  int n = x.size();
  Rcpp::NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    if(x[i] >= threshold) {
      out[i] = -std::pow(x[i], -2);
    } else {
      out[i] = -std::pow(threshold, -2);
    }
  }
  return out;
}

EL getEL(const arma::mat& g,
         const int& maxit,
         const double& abstol) {
  const int n = g.n_rows;
  const int p = g.n_cols;

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
    f0 = -Rcpp::sum(plog(Rcpp::wrap(arg), 1 / n));
    // J matrix & y vector
    arma::vec v1(Rcpp::sqrt(-d2plog(Rcpp::wrap(arg), 1 / n)));
    arma::vec v2(dplog(Rcpp::wrap(arg), 1 / n));
    J = g.each_col() % v1;
    y = v2 / v1;
    // update lambda by NR method with least square & step halving
    arma::qr_econ(Q, R, J);
    lc = l + arma::solve(arma::trimatu(R), Q.t() * y,
                         arma::solve_opts::fast);
    double alpha = 1;
    while(-Rcpp::sum(plog(Rcpp::wrap(1 + g * lc), 1 / n)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * arma::solve(arma::trimatu(R), Q.t() * y,
                                   arma::solve_opts::fast);
    }
    // update function value
    l = lc;
    arg = 1 + g * l;
    f1 = -Rcpp::sum(plog(Rcpp::wrap(arg), 1 / n));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      arma::vec v1(Rcpp::sqrt(-d2plog(Rcpp::wrap(arg), 1 / n)));
      arma::vec v2(dplog(Rcpp::wrap(arg), 1 / n));
      J = g.each_col() % v1;
      y = v2 / v1;
      convergence = true;
    } else {
      if(iterations == maxit) {
        break;
      }
      iterations++;
    }
  }

  EL result;
  result.nlogLR = -f1;
  result.lambda = l;
  result.gradient = -J.t() * y;
  result.iterations = iterations;
  result.convergence = convergence;
  return result;
}

double dfp1dcpp(arma::vec gr) {
  /* This function can be modified to for a general direction finding problem
     in higher dimensions
  */
  /* Prototype:
  arma::mat LHS = {{1, 0, 1}, {0, 1, -1}, {1, -1, 0}};
  arma::vec RHS = {0, 0, 0};
  RHS.subvec(0, 1) = -gr;
  arma::vec sol = solve(LHS, RHS);
  return sol(0);
  */
  return -sum(gr) / 2;
}

void theta2lambda_bibd(const arma::rowvec &theta,
                       arma::vec &lambda, bool &convex_hull_check,
                       const arma::mat &x, const arma::mat &c,
                       const int &n, const int &p,
                       int maxit, double abstol) {
  // estimating function
  arma::mat theta_mat = c.each_row() % theta;
  arma::mat g = x - theta_mat;

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
    f0 = -sum(plog(Rcpp::wrap(arg), 1 / n));
    // J matrix & y vector
    arma::vec v1(Rcpp::sqrt(-d2plog(Rcpp::wrap(arg), 1 / n)));
    arma::vec v2(dplog(Rcpp::wrap(arg), 1 / n));
    J = g.each_col() % v1;
    y = v2 / v1;
    // update lambda by NR method with least square & step halving
    arma::qr_econ(Q, R, J);
    lc = l + arma::solve(R, trans(Q) * y);
    double alpha = 1;
    while(-sum(plog(Rcpp::wrap(1 + g * lc), 1 / n)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * solve(R, trans(Q) * y);
    }
    // update function value
    l = lc;
    arg = 1 + g * l;
    f1 = -sum(plog(Rcpp::wrap(arg), 1 / n));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      arma::vec v1(Rcpp::sqrt(-d2plog(Rcpp::wrap(arg), 1 / n)));
      arma::vec v2(dplog(Rcpp::wrap(arg), 1 / n));
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
  lambda = l;
  convex_hull_check = convergence;
}

void lambda2theta_bibd(arma::rowvec &theta, arma::vec lambda,
                       const arma::mat &x, const arma::mat &c,
                       const std::vector<int> &pair, double gamma,
                       const int &n, const int &p,
                       int maxit, double abstol) {
  // estimating function
  arma::mat theta_mat = c.each_row() % theta;
  arma::mat g = x - theta_mat;

  // minimization
  arma::vec arg = 1 + g * lambda;
  arma::vec v(dplog(Rcpp::wrap(arg), 1 / n));
  arma::mat result2 = diagmat(v) * c;

  // gradient
  arma::rowvec gradient = -sum(result2, 0) % trans(lambda);

  // update parameter
  theta = theta - gamma * gradient;

  //projection
  double avg = (theta(pair[0] - 1) + theta((pair[1] - 1))) / 2;
  theta.at(pair[0] - 1) = avg;
  theta.at(pair[1] - 1) = avg;
}

std::vector<std::vector<int>> all_pairs(const int &p) {
  // initialize a vector of vectors
  std::vector<std::vector<int>> pairs;
  // the size of vector is p choose 2
  pairs.reserve(p * (p - 1) / 2);
  // fill in each elements(pairs)
  for(int i = 1; i < p + 1; i++)
  {
    for(int j = i + 1; j < p + 1; j++)
    {
      pairs.emplace_back(std::vector<int> {i, j});
    }
  }
  return pairs;
}

arma::mat cov_ibd(const arma::mat& x, const arma::mat& c) {
  // number of blocks
  int n = x.n_rows;
  // estimator(global minimizer)
  arma::vec theta_hat = n * arma::trans(arma::mean(x, 0) / arma::sum(c, 0));
  // estimating function
  arma::mat g = x - c.each_row() % theta_hat.t();
  // covariance estimate
  arma::mat vhat = (g.t() * (g)) / n;

  return vhat;
}

arma::mat g_mean(const arma::vec &theta,
                 arma::mat x) {
  // estimating function
  x.each_row() -= theta.t();

  return x;
}

arma::mat g_ibd(const arma::vec& theta,
                const arma::mat& x,
                const arma::mat& c) {
  return x - c.each_row() % theta.t();
}

arma::vec linear_projection(const arma::vec &theta,
                            const arma::mat &L,
                            const arma::vec &rhs) {
  return theta - L.t() * inv_sympd(L * L.t()) * (L * theta - rhs);
}

arma::vec lambda2theta_ibd(const arma::vec& lambda,
                           const arma::vec& theta,
                           const arma::mat& g,
                           const arma::mat& c,
                           const double& gamma) {
  arma::vec arg = 1 + g * lambda;
  arma::vec dplog_vec = dplog(Rcpp::wrap(arg), 1 / g.n_rows);
  // gradient
  arma::vec gradient = -arma::sum(arma::diagmat(dplog_vec) * c, 0).t() % lambda;
  // update theta by GD with lambda fixed
  arma::vec theta_hat = theta - gamma * gradient;

  return theta_hat;
}

double threshold_pairwise_ibd(const arma::mat& x,
                              const arma::mat& c,
                              const int& B,
                              const double& level) {
  /// parameters ///
  const int p = x.n_cols;   // number of points(treatments)
  const int m = p * (p - 1) / 2;    // number of hypotheses
  const arma::vec conf_level = {level};    // confidence level

  /// A hat matrices ///
  arma::cube A_hat(p, p, m);
  // covariance estimate
  const arma::mat V_hat = cov_ibd(x, c);
  // vector of pairs
  const std::vector<std::vector<int>> pairs = all_pairs(p);
  for (int i = 0; i < m; i++) {
    arma::rowvec R = arma::zeros(1, p);
    R(pairs[i][0] - 1) = 1;
    R(pairs[i][1] - 1) = -1;
    A_hat.slice(i) = (R.t() * R) / arma::as_scalar(R * V_hat * R.t());
  }

  // U hat matrices
  const arma::mat U_hat = arma::mvnrnd(arma::zeros(p), V_hat, B);

  // B bootstrap replicates(B x m matrix)
  arma::mat bootstrap_sample(B, m);
  for (int i = 0; i < m; i++) {
    bootstrap_sample.col(i) =
      arma::diagvec(arma::trans(U_hat) * A_hat.slice(i) * U_hat);
  }

  return
    arma::as_scalar(arma::quantile(arma::max(bootstrap_sample, 1), conf_level));
}

