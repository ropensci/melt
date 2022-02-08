#include "EL.h"

// [[Rcpp::export]]
Rcpp::List EL_mean(const Eigen::Map<Eigen::VectorXd>& par,
                   const Eigen::Map<Eigen::MatrixXd>& x,
                   const int maxit,
                   const double abstol) {
  // check 'par' and 'x'
  const int n = x.rows();
  const int p = x.cols();
  if (par.size() != p) {
    Rcpp::stop("dimensions of 'par' and 'x' do not match");
  }
  if (n < 2) {
    Rcpp::stop("not enough 'x' observations");
  }
  const Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(x);
  if (lu_decomp.rank() != p || n <= p) {
    Rcpp::stop("matrix 'x' must have full column rank");
  }
  const Eigen::VectorXd estimate = x.colwise().mean();
  const EL2 el(par, x, "mean", maxit, abstol);
  // const EL el(x.rowwise() - par.transpose(), p * 100, maxit, abstol);

  const double chisq_statistic = 2 * el.nlogLR;
  Rcpp::Function pchisq("pchisq");
  const double pval = Rcpp::as<double>(
    pchisq(chisq_statistic, Rcpp::Named("df") = p,
           Rcpp::Named("lower.tail") = false));
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("optim") = Rcpp::List::create(
      Rcpp::Named("type") = "mean",
      Rcpp::Named("lambda") = el.lambda,
      Rcpp::Named("logLR") = -el.nlogLR,
      Rcpp::Named("iterations") = el.iterations,
      Rcpp::Named("convergence") = el.convergence),
    Rcpp::Named("statistic") = 2 * el.nlogLR,
    Rcpp::Named("df") = p,
    Rcpp::Named("p.value") = pval,
    Rcpp::Named("estimate") = estimate,
    Rcpp::Named("null.value") = par,
    Rcpp::Named("alternative") = "two.sided",
    Rcpp::Named("method") = "One sample EL test");
  result.attr("class") = Rcpp::CharacterVector({"el_test"});
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector EL_confidence_interval_mean(
    const Eigen::Map<Eigen::MatrixXd>& x,
    const double init,
    const double cutoff,
    const int maxit = 100,
    const double abstol = 1e-8) {
  // upper endpoint
  double upper_lb = init;
  double upper_size = 1;
  double upper_ub = init + upper_size;
  // upper bound for upper endpoint
  while (2 * EL(x.array() - upper_ub, 100, maxit, abstol).nlogLR <= cutoff) {
    upper_lb = upper_ub;
    upper_ub += upper_size;
  }
  // approximate upper bound by numerical search
  while (upper_ub - upper_lb > 1e-04) {
    if (2 * EL(x.array() - ((upper_lb + upper_ub) / 2), 100,
               maxit, abstol).nlogLR > cutoff) {
      upper_ub = (upper_lb + upper_ub) / 2;
    } else {
      upper_lb = (upper_lb + upper_ub) / 2;
    }
  }

  // lower endpoint
  double lower_ub = init;
  double lower_size = 1;
  double lower_lb = init - lower_size;


  // lower bound for lower endpoint
  while (2 * EL(x.array() - lower_lb, 100, maxit, abstol).nlogLR <= cutoff) {
    lower_ub = lower_lb;
    lower_lb -= lower_size;
  }
  // approximate lower bound by numerical search
  while (lower_ub - lower_lb > 1e-04) {
    if (2 * EL(x.array() - ((lower_lb + lower_ub) / 2), 100,
               maxit, abstol).nlogLR > cutoff) {
      lower_lb = (lower_lb + lower_ub) / 2;
    } else {
      lower_ub = (lower_lb + lower_ub) / 2;
    }
  }
  Rcpp::NumericVector v = {lower_ub, upper_lb};
  return v;
}

// [[Rcpp::export]]
Eigen::MatrixXd EL_confidence_interval(
    const Eigen::Map<Eigen::MatrixXd>& x,
    const std::string type,
    const Eigen::Map<Eigen::VectorXd>& init,
    const double cutoff,
    const int maxit = 100,
    const double abstol = 1e-8) {
  const int p = x.cols();
  Eigen::MatrixXd ci(p, 2);
  for (int j = 0; j < p; ++j) {
    // upper endpoint
    double upper_lb = init[j];
    double upper_size = 1;
    double upper_ub = init[j] + upper_size;
    // upper bound for upper endpoint
    while (2 * EL2(Eigen::Matrix<double, 1, 1>(upper_ub), x, type, maxit,
                   abstol).nlogLR <= cutoff) {
      upper_lb = upper_ub;
      upper_ub += upper_size;
    }
    // approximate upper bound by numerical search
    while (upper_ub - upper_lb > 1e-04) {
      if (2 * EL2(Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2), x,
                  type, maxit, abstol).nlogLR > cutoff) {
        upper_ub = (upper_lb + upper_ub) / 2;
      } else {
        upper_lb = (upper_lb + upper_ub) / 2;
      }
    }
    // lower endpoint
    double lower_ub = init[j];
    double lower_size = 1;
    double lower_lb = init[j] - lower_size;

    // lower bound for lower endpoint
    while (2 * EL2(Eigen::Matrix<double, 1, 1>(lower_lb), x, type, maxit,
                   abstol).nlogLR <= cutoff) {
      lower_ub = lower_lb;
      lower_lb -= lower_size;
    }
    // approximate lower bound by numerical search
    while (lower_ub - lower_lb > 1e-04) {
      if (2 * EL2(Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2), x,
                  type, maxit, abstol).nlogLR > cutoff) {
        lower_lb = (lower_lb + lower_ub) / 2;
      } else {
        lower_ub = (lower_lb + lower_ub) / 2;
      }
    }
    ci(j, 0) = lower_lb;
    ci(j, 1) = upper_ub;
  }
  // // upper endpoint
  // double upper_lb = init[0];
  // double upper_size = 1;
  // double upper_ub = init[0] + upper_size;
  // // upper bound for upper endpoint
  // while (2 * EL2(Eigen::Matrix<double, 1, 1>(upper_ub), x, type, maxit,
  //                abstol).nlogLR <= cutoff) {
  //   upper_lb = upper_ub;
  //   upper_ub += upper_size;
  // }
  // // approximate upper bound by numerical search
  // while (upper_ub - upper_lb > 1e-04) {
  //   if (2 * EL(x.array() - ((upper_lb + upper_ub) / 2), 100,
  //              maxit, abstol).nlogLR > cutoff) {
  //     upper_ub = (upper_lb + upper_ub) / 2;
  //   } else {
  //     upper_lb = (upper_lb + upper_ub) / 2;
  //   }
  // }
  //
  // // lower endpoint
  // double lower_ub = init[0];
  // double lower_size = 1;
  // double lower_lb = init[0] - lower_size;
  //
  // // lower bound for lower endpoint
  // while (2 * EL(x.array() - lower_lb, 100, maxit, abstol).nlogLR <= cutoff) {
  //   lower_ub = lower_lb;
  //   lower_lb -= lower_size;
  // }
  // // approximate lower bound by numerical search
  // while (lower_ub - lower_lb > 1e-04) {
  //   if (2 * EL(x.array() - ((lower_lb + lower_ub) / 2), 100,
  //              maxit, abstol).nlogLR > cutoff) {
  //     lower_lb = (lower_lb + lower_ub) / 2;
  //   } else {
  //     lower_ub = (lower_lb + lower_ub) / 2;
  //   }
  // }
  // Rcpp::NumericVector v = {lower_ub, upper_lb};
  return ci;
}

