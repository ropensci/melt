#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]
arma::vec plog(const arma::vec& x, double threshold) {
  const int n = x.n_elem;
  arma::vec out(n);
  for(int i = 0; i < n; ++i) {
    if(x(i) >= threshold) {
      out(i) = std::log(x(i));
    } else {
      out(i) = std::log(threshold) - 1.5 + 2 * std::pow(threshold, -1) * x(i) -
        std::pow(x(i) / threshold, 2) / 2;
    }
  }
  return out;
}
Rcpp::NumericVector plog_old(Rcpp::NumericVector x, double threshold) {
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

arma::vec dplog(const arma::vec& x, double threshold) {
  const int n = x.n_elem;
  arma::vec out(n);
  for(int i = 0; i < n; ++i) {
    if(x(i) >= threshold) {
      out(i) = std::pow(x(i), -1);
    } else {
      out(i) = 2 * std::pow(threshold, -1) - x(i) * std::pow(threshold, -2);
    }
  }
  return out;
}
Rcpp::NumericVector dplog_old(Rcpp::NumericVector x, double threshold) {
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

arma::vec d2plog(const arma::vec& x, double threshold) {
  const int n = x.n_elem;
  arma::vec out(n);
  for(int i = 0; i < n; ++i) {
    if(x(i) >= threshold) {
      out(i) = -std::pow(x(i), -2);
    } else {
      out(i) = -std::pow(threshold, -2);
    }
  }
  return out;
}
Rcpp::NumericVector d2plog_old(Rcpp::NumericVector x, double threshold) {
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
    f0 = -arma::sum(plog(arg, 1.0 / n));
    // J matrix & y vector
    arma::vec v1 = arma::sqrt(-d2plog(arg, 1.0 / n));
    arma::vec v2 = dplog(arg, 1.0 / n);
    J = g.each_col() % v1;
    y = v2 / v1;
    // update lambda by NR method with least square & step halving
    arma::qr_econ(Q, R, J);
    lc = l + arma::solve(arma::trimatu(R), Q.t() * y,
                         arma::solve_opts::fast);
    double alpha = 1.0;
    while(-arma::sum(plog(1 + g * lc, 1.0 / n)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * arma::solve(arma::trimatu(R), Q.t() * y,
                                   arma::solve_opts::fast);
    }
    // update function value
    l = lc;
    arg = 1 + g * l;
    f1 = -arma::sum(plog(arg, 1.0 / n));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      arma::vec v1 = arma::sqrt(-d2plog(arg, 1.0 / n));
      arma::vec v2 = dplog(arg, 1.0 / n);
      J = g.each_col() % v1;
      y = v2 / v1;
      convergence = true;
    } else {
      if(iterations == maxit) {
        break;
      }
      ++iterations;
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

std::vector<std::array<int, 2>> all_pairs(const int p) {
  // initialize a vector of vectors
  std::vector<std::array<int, 2>> pairs;
  // the size of vector is p choose 2
  pairs.reserve(p * (p - 1) / 2);
  // fill in each elements(pairs)
  for(int i = 1; i < p + 1; ++i)
  {
    for(int j = i + 1; j < p + 1; ++j)
    {
      pairs.emplace_back(std::array<int, 2>{j, i});
      // pairs.emplace_back(std::vector<int> {i, j});
    }
  }
  return pairs;
}

arma::mat g_mean(const arma::vec &theta,
                 arma::mat x) {
  // estimating function
  x.each_row() -= theta.t();

  return x;
}

arma::vec linear_projection(const arma::vec &theta,
                            const arma::mat &L,
                            const arma::vec &rhs) {
  return theta - L.t() * inv_sympd(L * L.t()) * (L * theta - rhs);
}

arma::mat bootstrap_sample(const arma::mat& x)
{
  const int n = x.n_rows;
  const int p = x.n_cols;
  // uniform random sample of index integers
  arma::vec boostrap_index =
    arma::randi<arma::vec>(n, arma::distr_param(0, n - 1));

  // resampling with replacement
  arma::mat bootstrap_sample(n, p);
  for(int i = 0; i < n; ++i) {
    bootstrap_sample.row(i) = x.row(boostrap_index(i));
  }

  return bootstrap_sample;
}
