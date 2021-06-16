#include "utils.h"

arma::vec plog(const arma::vec& x) {
  const unsigned int n = x.n_elem;
  const double a1 = -std::log(n) - 1.5;
  const double a2 = 2.0 * n;
  const double a3 = 0.5 * n * n;
  arma::vec out(n);
  for (unsigned int i = 0; i < n; ++i) {
    if (n * x.at(i) < 1) {
      out.at(i) = a1 + a2 * x.at(i) - a3 * x.at(i) * x.at(i);
    } else {
      out.at(i) = std::log(x.at(i));
    }
  }
  return out;
}
arma::vec plog_old(const arma::vec& x, const double threshold) {
  arma::vec out(x.n_elem);
  for (unsigned int i = 0; i < x.n_elem; ++i) {
    if (x.at(i) < threshold) {
      out.at(i) = std::log(threshold) - 1.5 + 2 * std::pow(threshold, -1) * x.at(i) -
        std::pow(x.at(i) / threshold, 2) / 2;
    } else {
      out.at(i) = std::log(x.at(i));
    }
  }
  return out;
}

arma::vec dplog(const arma::vec& x) {
  const unsigned int n = x.n_elem;
  const double a1 = 2.0 * n;
  const double a2 = 1.0 * n * n;
  arma::vec out(n);
  for (unsigned int i = 0; i < n; ++i) {
    if (n * x.at(i) < 1) {
      out.at(i) = a1 - a2 * x.at(i);
    } else {
      out.at(i) = 1.0 / x.at(i);
    }
  }
  return out;
}
arma::vec dplog_old(const arma::vec& x, const double threshold) {
  arma::vec out(x.n_elem);
  for (unsigned int i = 0; i < x.n_elem; ++i) {
    if (x.at(i) < threshold) {
      out.at(i) = 2 * std::pow(threshold, -1) - x.at(i) * std::pow(threshold, -2);
    } else {
      out.at(i) = std::pow(x.at(i), -1);
    }
  }
  return out;
}

arma::vec d2plog(const arma::vec& x, const double threshold) {
  arma::vec out(x.n_elem);
  for (unsigned int i = 0; i < x.n_elem; ++i) {
    if (x.at(i) < threshold) {
      out.at(i) = -std::pow(threshold, -2);
    } else {
      out.at(i) = -std::pow(x.at(i), -2);
    }
  }
  return out;
}
arma::vec d2plog_old(const arma::vec& x, const double threshold) {
  const int n = x.n_elem;
  arma::vec out(n);
  for (int i = 0; i < n; ++i) {
    if (x(i) >= threshold) {
      out(i) = -std::pow(x(i), -2);
    } else {
      out(i) = -std::pow(threshold, -2);
    }
  }
  return out;
}


EL getEL(const arma::mat& g,
         const int maxit,
         const double abstol) {
  const int n = g.n_rows;
  const int p = g.n_cols;

  // minimization
  arma::vec l(p, arma::fill::zeros);
  arma::vec lc(p);
  arma::vec arg = 1 + g * l;
  arma::vec y(n);
  arma::mat J(n, p);
  arma::mat Q(n, n);
  arma::mat R(p, p);
  double f0;
  double f1;
  int iterations = 0;
  bool convergence = false;
  while (!convergence) {
    // function evaluation(initial)
    f0 = -arma::sum(plog(arg));
    // J matrix & y vector
    arma::vec v1 = arma::sqrt(-d2plog(arg, 1.0 / n));
    arma::vec v2 = dplog(arg);
    J = g.each_col() % v1;
    y = v2 / v1;
    // update lambda by NR method with least square & step halving
    arma::qr_econ(Q, R, J);
    arma::vec update =
      arma::solve(arma::trimatu(R), Q.t() * y, arma::solve_opts::fast);
    lc = l + update;
    double alpha = 1.0;
    while (-arma::sum(plog(1 + g * lc)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * update;
    }
    // update function value
    l = lc;
    arg = 1 + g * l;
    f1 = -arma::sum(plog(arg));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      arma::vec v1 = arma::sqrt(-d2plog(arg, 1.0 / n));
      arma::vec v2 = dplog(arg);
      J = g.each_col() % v1;
      y = v2 / v1;
      convergence = true;
    } else {
      if (iterations == maxit) {
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
  for (int i = 1; i < p + 1; ++i) {
    for (int j = i + 1; j < p + 1; ++j) {
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

arma::vec linear_projection(const arma::vec& theta,
                            const arma::mat& L,
                            const arma::vec& rhs) {
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
  for (int i = 0; i < n; ++i) {
    bootstrap_sample.row(i) = x.row(boostrap_index(i));
  }

  return bootstrap_sample;
}
