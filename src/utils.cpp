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
                       Rcpp::IntegerVector pair, double gamma,
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
  double avg = (theta(pair(0) - 1) + theta((pair(1) - 1))) / 2;
  theta.at(pair(0) - 1) = avg;
  theta.at(pair(1) - 1) = avg;
}

