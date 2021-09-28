#include "EL.h"

EL2::EL2(const Eigen::Ref<const Eigen::MatrixXd>& g,
         const double threshold,
         const int maxit,
         const double abstol) {
  // maximization
  lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  iterations = 0;
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

// EL getEL(const Eigen::Ref<const Eigen::MatrixXd>& g,
//          const int maxit,
//          const double abstol) {
//   // maximization
//   // Eigen::VectorXd lambda = Eigen::VectorXd::Zero(g.cols());
//   Eigen::VectorXd lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
//   double f1;
//   int iterations = 0;
//   bool convergence = false;
//   while (!convergence && iterations != maxit) {
//     // plog class
//     PSEUDO_LOG log_tmp(Eigen::VectorXd::Ones(g.rows()) + g * lambda);
//
//     // // function evaluation
//     // const double f0 = log_tmp.plog_sum;
//
//     // J matrix
//     const Eigen::MatrixXd J = g.array().colwise() * log_tmp.sqrt_neg_d2plog;
//
//     // // J matrix & y vector
//     // Eigen::ArrayXd v1 = log_tmp.sqrt_neg_d2plog;
//     // Eigen::ArrayXd v2 = log_tmp.dplog;
//     // Eigen::MatrixXd J = g.array().colwise() * v1;
//     /// Eigen::VectorXd y = v2 / v1;
//
//     // prpose new lambda by NR method with least square
//     Eigen::VectorXd&& step =
//       (J.transpose() * J).ldlt().solve(
//           J.transpose() * (log_tmp.dplog / log_tmp.sqrt_neg_d2plog).matrix());
//
//     // update function value
//     f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
//
//     // step halving to ensure validity
//     while (f1 < log_tmp.plog_sum) {
//       step /= 2;
//       f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * (lambda + step));
//     }
//
//     // update lambda
//     lambda += step;
//
//     // convergence check
//     if (f1 - log_tmp.plog_sum < abstol) {
//       // Eigen::ArrayXd v1 = log_tmp.sqrt_neg_d2plog;
//       // Eigen::ArrayXd v2 = log_tmp.dplog;
//       // Eigen::MatrixXd J = g.array().colwise() * v1;
//       // Eigen::VectorXd y = v2 / v1;
//       convergence = true;
//     } else {
//       ++iterations;
//     }
//   }
//
//   return {lambda, f1, iterations, convergence};
// }
