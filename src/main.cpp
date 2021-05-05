#include <RcppArmadillo.h>
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
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
List el_mean(arma::rowvec theta, arma::mat x,
               int maxit = 100, double abstol = 1e-8) {
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
    f0 = -sum(plog(wrap(arg), 1 / n));
    // J matrix & y vector
    arma::vec v1(Rcpp::sqrt(-d2plog(wrap(arg), 1 / n)));
    arma::vec v2(dplog(wrap(arg), 1 / n));
    J = g.each_col() % v1;
    y = v2 / v1;
    // update lambda by NR method with least square & step halving
    arma::qr_econ(Q, R, J);
    lc = l + arma::solve(R, trans(Q) * y);
    double alpha = 1;
    while(-sum(plog(wrap(1 + g * lc), 1 / n)) > f0) {
      alpha = alpha / 2;
      lc = l + alpha * solve(R, trans(Q) * y);
    }
    // update function value
    l = lc;
    arg = 1 + g * l;
    f1 = -sum(plog(wrap(arg), 1 / n));
    // convergence check & parameter update
    if (f0 - f1 < abstol) {
      arma::vec v1(Rcpp::sqrt(-d2plog(wrap(arg), 1 / n)));
      arma::vec v2(dplog(wrap(arg), 1 / n));
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
    Rcpp::Named("lambda") = NumericVector(l.begin(), l.end()),
    Rcpp::Named("grad") = as<std::vector<double>>(wrap(-trans(J) * y)),
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("convergence") = convergence);
}


//' Two sample test for equal mean
//'
//' Two sample test for equal mean
//'
//' @param x a vector of data for one  group.
//' @param y a vector of data for the other  group.
//' @param b a momentum parameter for minimization. Defaults to .1.
//' @param alpha an optional step size. Defaults to 1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @examples
//' x <- rnorm(100)
//' y <- rnorm(100)
//' test2sample2_cpp(x, y)
//'
//' @export
// [[Rcpp::export]]
List test2sample2_cpp(NumericVector x, NumericVector y, double b = .9, double alpha = 1,
                      unsigned int maxit = 1000, double abstol = 1e-8) {
  List result;
  vec sample_mean = {mean(x), mean(y)};
  vec ub_vec = {max(x), max(y), max(sample_mean)};
  vec lb_vec = {min(x), min(y), min(sample_mean)};
  double ub = min(ub_vec);
  double lb = max(lb_vec);
  if(ub <= lb) {
    result["nlogLR"] = datum::inf;;
    result["convergence"] = -1;
    return result;
  }

  // initialization
  double par = (lb + ub) / 2;
  unsigned int nx = x.size();
  unsigned int ny = y.size();
  unsigned int N = std::max(nx, ny);
  unsigned int iterations{};
  int convergence{};
  double v{};
  double lx;
  double ly;
  arma::vec grad;
  double d;


  // minimization
  while (convergence == 0) {
    // lambda update
    lx = el_mean({par}, x)["lambda"];
    ly = el_mean({par}, y)["lambda"];
    // gradient
    grad = {sum(dplog(1 + lx * (x - par), 1 / nx)) * (-lx) / N,
            sum(dplog(1 + ly * (y - par), 1 / ny)) * (-ly) / N};
    // direction
    d = dfp1dcpp(grad);
    // direction change reverts momentum
    if (sign(d) != sign(v)) {
      v = 0; alpha = alpha / 2;
    }
    // lb, ub update
    if (sign(d) > 0) {
      lb = par;
    } else {
      ub = par;
    }
    // convergence check & parameter update
    if ((std::abs(d * sum(grad)) < abstol || ub - lb < abstol) && iterations > 0) {
      convergence = 1;
    } else {
      iterations++;
      // step halving to satisfy convex hull constraint
      v = b * v + d;
      while (par + alpha * v <= lb || par + alpha * v >= ub) {
        alpha = alpha / 2;
      }
      par = par + alpha * v;
      if(iterations == maxit) {
        break;
      }
    }
  }

  // result
  return Rcpp::List::create(
    Rcpp::Named("par") = par,
    Rcpp::Named("nlogLR") = sum(Rcpp::log(1 + (x - par) * lx)) +
      sum(Rcpp::log(1 + (y - par) * ly)),
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("convergence") = convergence
  );
}


//' Two sample test for equal mean
//'
//' Two sample test for equal mean
//'
//' @param x a vector of data for one  group.
//' @param y a vector of data for the other  group.
//' @param b a momentum parameter for minimization. Defaults to .1.
//' @param alpha an optional step size. Defaults to 1.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @examples
//' x <- rnorm(100)
//' y <- rnorm(100)
//' test2sample777_cpp(x, y)
//'
//' @export
// [[Rcpp::export]]
List test2sample777_cpp(NumericVector x, NumericVector y, double b = .9, double alpha = 1,
                        unsigned int maxit = 1000, double abstol = 1e-8) {
  List result;
  vec sample_mean = {mean(x), mean(y)};
  vec ub_vec = {max(x), max(y), max(sample_mean)};
  vec lb_vec = {min(x), min(y), min(sample_mean)};
  double ub = min(ub_vec);
  double lb = max(lb_vec);
  if(ub <= lb) {
    result["nlogLR"] = datum::inf;;
    result["convergence"] = -1;
    return result;
  }

  // initialization
  double par = (lb + ub) / 2;
  unsigned int nx = x.size();
  unsigned int ny = y.size();
  unsigned int N = std::max(nx, ny);
  unsigned int iterations{};
  int convergence{};
  double v{};
  double lx;
  double ly;
  arma::vec grad;
  double d;

  // linearization
  double par0 = par;
  double lx0 = el_mean({par}, x)["lambda"];
  double ly0 = el_mean({par}, y)["lambda"];
  double mx = -sum(1 / pow(1 + lx0 * (x - par), 2)) /
    sum(pow(x - par, 2) / pow(1 + lx0 * (x - par), 2));
  double my = -sum(1 / pow(1 + ly0 * (y - par) , 2)) /
    sum(pow(y - par ,2) / pow(1 + ly0 * (y - par) , 2));

  // minimization
  while (convergence == 0) {
    // lambda update
    lx = mx * (par - par0) + lx0;
    ly = my * (par - par0) + ly0;
    // gradient
    grad = {sum(dplog(1 + lx * (x - par), 1 / nx)) * (-lx) / N,
            sum(dplog(1 + ly * (y - par), 1 / ny)) * (-ly) / N};
    // direction
    d = dfp1dcpp(grad);
    // direction change reverts momentum
    if (sign(d) != sign(v)) {
      v = 0; alpha = alpha / 2;
    }
    // lb, ub update
    if (sign(d) > 0) {
      lb = par;
    } else {
      ub = par;
    }
    // convergence check & parameter update
    if ((std::abs(d * sum(grad)) < abstol || ub - lb < abstol) && iterations > 0) {
      convergence = 1;
    } else {
      iterations++;
      // step halving to satisfy convex hull constraint
      v = b * v + d;
      while (par + alpha * v <= lb || par + alpha * v >= ub) {
        alpha = alpha / 2;
      }
      par = par + alpha * v;
      if(iterations == maxit) {
        break;
      }
    }
  }

  // result
  return Rcpp::List::create(
    Rcpp::Named("par") = par,
    Rcpp::Named("nlogLR") = sum(Rcpp::log(1 + (x - par) * lx)) +
      sum(Rcpp::log(1 + (y - par) * ly)),
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("convergence") = convergence
  );
}


//' Two sample test for BIBD
//'
//' Two sample test for BIBD
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param pair a pair of index for parameters to be tested.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @export
// [[Rcpp::export]]
List test_pair(const arma::mat &x,
               const arma::mat &c,
               const std::vector<int> &pair,
               int maxit = 1000,
               double abstol = 1e-8) {
  /// index pair validation ///
  int p = x.n_cols;   // number of points(parameters)
  // vector for entire points
  Rcpp::IntegerVector points = Rcpp::seq(1, p);
  // union of pair and points
  Rcpp::IntegerVector pair_Rcpp = Rcpp::wrap(pair);
  Rcpp::IntegerVector total = Rcpp::union_(pair_Rcpp, points);
  if (pair.size() != 2 ||
      sum(duplicated(pair_Rcpp)) ||
      !Rcpp::setequal(points, total)) {
      stop("'pair' must be a vector of two distinct integers ranging from one to the total number of points.");
  }

  /// convex hull check(necessary condition) ///
  // upper bound for the common parameter value
  double ub = std::min(arma::max(x.col(pair[0] - 1)),
                       arma::max(x.col(pair[1] - 1)));
  // lower bound form the commom parameter value
  double lb = std::max(arma::min(x.col(pair[0] - 1)),
                       arma::min(x.col(pair[1] - 1)));
  if (ub <= lb) {
    return Rcpp::List::create(
      Rcpp::Named("nlogLR") = datum::inf,
      Rcpp::Named("convergence") = -1);
  }

  /// initialization ///
  int n = x.n_rows;
  // initial parameter value set as column means multiplied by b / r
  arma::rowvec theta = n / accu(x.col(0) != 0) * mean(x, 0);
  // equality constraint imposed on the initial value
  double avg = (theta(pair[0] - 1) + theta((pair[1] - 1))) / 2;
  theta.at(pair[0] - 1) = avg;
  theta.at(pair[1] - 1) = avg;
  bool convex_hull_check;   // for checking convex hull constraint
  double f0;    // for current function value(-logLR)
  double f1;    // for updated function value
  arma::vec arg;    // for expression evaluated insided the pseudo log function
  arma::vec lambda;   // initial lambda
  theta2lambda_bibd(theta, lambda, convex_hull_check, x, c, n, p, 100, 1e-8);
  // If the convex hull constraint is not satisfied at the initial value, stop.
  // Otherwise, compute initial function value.
  if (!convex_hull_check) {
    stop("convex hull constraint not satisfied at the initial value.");
  } else {
    arg = 1 + (x - c.each_row() % theta) * lambda;
    f0 = sum(plog(wrap(arg), 1 / n));
    f1 = f0;
  }

  /// minimization(projected gradient descent) ///
  double gamma = 2 * pow(n, -1);    // step size
  int convergence{};
  int iterations = 0;
  while (convergence == 0) {
    if (f0 - f1 < abstol && iterations > 0) {
      convergence = 1;
    } else {
      // update parameter
      lambda2theta_bibd(theta, lambda, x, c, pair, gamma, n, p, 100, 1e-8);
      // update lambda
      theta2lambda_bibd(theta, lambda, convex_hull_check, x, c, n, p, 100, 1e-8);
      // update function value
      f0 = f1;
      arg = 1 + (x - c.each_row() % theta) * lambda;
      f1 = sum(plog(wrap(arg), 1 / n));
      // step halving to ensure that the updated function value be
      // strinctly less than the current function value
      while (f0 <= f1) {
        gamma = gamma / 2;
        lambda2theta_bibd(theta, lambda, x, c, pair, gamma, n, p, 100, 1e-8);
        arg = 1 + (x - c.each_row() % theta) * lambda;
        f1 = sum(plog(wrap(arg), 1 / n));
      }
      iterations++;
      if (iterations == maxit) {
        break;
      }
    }
  }

  // result
  return Rcpp::List::create(
    Rcpp::Named("theta") = NumericVector(theta.begin(), theta.end()),
    Rcpp::Named("lambda") = NumericVector(lambda.begin(), lambda.end()),
    Rcpp::Named("nlogLR") = f1,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("convergence") = convergence);
}

//' Pairwise comparison for BIBD
//'
//' Pairwise comparison for BIBD
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
//' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
//'
//' @export
// [[Rcpp::export]]
List pairwise_test(const arma::mat &x,
                   const arma::mat &c,
                   int maxit = 1000,
                   double abstol = 1e-8) {
  // number of points(parameters)
  int p = x.n_cols;
  // all pairs
  std::vector<std::vector<int>> pairs = all_pairs(p);
  // number of hypotheses
  int m = pairs.size();
  // test statistics(-2logLR)
  Rcpp::NumericVector pairs_2nlogLR(m);
  for (unsigned i = 0; i < m; i++) {
    double nlogLR = test_pair(x, c, pairs[i], maxit, abstol)["nlogLR"];
    pairs_2nlogLR(i) = 2 * nlogLR;
  }
  // std::vector<double> pairs_2nlogLR(m);
  // for (int i = 0; i < m; i++) {
  //   double nlogLR = test_pair(x, c, pairs[i], maxit, abstol)["nlogLR"];
  //   pairs_2nlogLR[i] = 2 * nlogLR;
  // }
  return Rcpp::List::create(
    Rcpp::Named("statistics") = pairs_2nlogLR);
}

//' Computes threshold of pairwise comparison for BIBD
//'
//' Computes threshold using parametric bootstrap
//'
//' @param x a matrix of data .
//' @param c an incidence matrix.
//' @param B an integer value for the number of bootstrap replicates.
//' @param alpha a significance level.
//'
//' @export
// [[Rcpp::export]]
List pairwise_threshold(const arma::mat &x,
                        const arma::mat &c,
                        const int &B,
                        const double &alpha) {
  /// parameters ///
  const int n = x.n_rows;   // number of blocks
  const int p = x.n_cols;   // number of points(treatments)
  const int m = p * (p - 1) / 2;    // number of hypotheses
  const arma::vec level = {1 - alpha};    // confidence level

  /// A hat matrices ///
  arma::cube A_hat(p, p, m);
  // covariance estimate
  const arma::mat V_hat = cov_estimator(x, c);
  // vector of pairs
  const std::vector<std::vector<int>> pairs = all_pairs(p);
  for (unsigned i = 0; i < m; i++) {
    arma::rowvec R(p);
    R.fill(0);
    R(pairs[i][0] - 1) = 1;
    R(pairs[i][1] - 1) = -1;
    A_hat.slice(i) = (trans(R) * R) / as_scalar(R * V_hat * trans(R));
  }

  // U hat matrices
  const arma::mat U_hat = arma::mvnrnd(arma::zeros(p), V_hat, B);

  // B bootstrap replicates(B x m matrix)
  arma::mat bootstrap_sample(B, m);
  for (unsigned i = 0; i < m; i++) {
    bootstrap_sample.col(i) =
      arma::diagvec(arma::trans(U_hat) * A_hat.slice(i) * U_hat);
  }

  return Rcpp::List::create(
    Rcpp::Named("threshold") =
      arma::as_scalar(arma::quantile(arma::max(bootstrap_sample, 1), level)));
}
