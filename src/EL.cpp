#include "EL.h"

EL::EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
       const double threshold,
       const int maxit,
       const double abstol) {
  // maximization
  lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  iterations = 1;
  convergence = false;
  while (!convergence && iterations != maxit) {
    // plog class
    PSEUDO_LOG log_tmp(Eigen::VectorXd::Ones(g.rows()) + g * lambda);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * log_tmp.sqrt_neg_d2plog;
    // prpose new lambda by NR method with least square
    Eigen::VectorXd step =
      (J.transpose() * J).ldlt().solve(
          J.transpose() * (log_tmp.dplog / log_tmp.sqrt_neg_d2plog).matrix());
    // update function value
    nlogLR =
      PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
    // step halving to ensure validity
    while (nlogLR < log_tmp.plog_sum) {
      step /= 2;
      nlogLR =
        PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
    }
    // update lambda
    lambda += step;

    // check convex hull constraint(stop if larger than threshold)
    if (nlogLR > threshold) {
      break;
    }

    // convergence check
    if (nlogLR - log_tmp.plog_sum < abstol) {
      convergence = true;
    } else {
      ++iterations;
    }
  }
}

EL2::EL2(const Eigen::Ref<const Eigen::VectorXd>& par,
         const Eigen::Ref<const Eigen::MatrixXd>& x,
         const std::string type,
         const int maxit,
         const double abstol) {
  std::map<std::string,
           std::function<Eigen::MatrixXd(
             const Eigen::Ref<const Eigen::VectorXd>&,
             const Eigen::Ref<const Eigen::MatrixXd>&)>> funcMap{
               {{"mean", g_mean},
               {"lm", g_mean}}
             };
  // Eigen::MatrixXd (*fcnPtr)(
  //     const Eigen::Ref<const Eigen::VectorXd>& par,
  //     const Eigen::Ref<const Eigen::MatrixXd>& x){ &g_mean };

  Eigen::MatrixXd g = funcMap[type](par, x);
  // maximization
  lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  iterations = 1;
  convergence = false;
  while (!convergence && iterations != maxit) {
    // plog class
    PSEUDO_LOG log_tmp(Eigen::VectorXd::Ones(g.rows()) + g * lambda);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * log_tmp.sqrt_neg_d2plog;
    // prpose new lambda by NR method with least square
    Eigen::VectorXd step =
      (J.transpose() * J).ldlt().solve(
          J.transpose() * (log_tmp.dplog / log_tmp.sqrt_neg_d2plog).matrix());
    // update function value
    nlogLR =
      PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
    // step halving to ensure validity
    while (nlogLR < log_tmp.plog_sum) {
      step /= 2;
      nlogLR =
        PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
    }
    // update lambda
    lambda += step;

    // // check convex hull constraint(stop if larger than threshold)
    // if (nlogLR > threshold) {
    //   break;
    // }

    // convergence check
    if (nlogLR - log_tmp.plog_sum < abstol) {
      convergence = true;
    } else {
      ++iterations;
    }
  }
}
