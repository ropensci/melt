// #include "mht-utils.h"
//
// Rcpp::List mht_(const std::string method,
//                 const Eigen::Map<Eigen::VectorXd>& est,
//                 const Eigen::Map<Eigen::MatrixXd>& x,
//                 const Eigen::Map<Eigen::VectorXd>& rhs,
//                 const Eigen::Map<Eigen::MatrixXd>& lhs,
//                 const int maxit,
//                 const int maxit_l,
//                 const double tol,
//                 const double tol_l,
//                 const Rcpp::Nullable<double> step,
//                 const Rcpp::Nullable<double> th,
//                 const int nthreads,
//                 const Eigen::Map<Eigen::ArrayXd>& wt,
//                 const Eigen::Map<Eigen::VectorXi>& q,
//                 const int m,
//                 const double level,
//                 const int B)
// {
//   // test statistics
//
//   // critical value
//   Rcpp::List result = Rcpp::List::create(
//     Rcpp::Named("cv") = cv_mvchisq(method, est, x, lhs, q, m, level, B, nthreads));
//   return result;
// }
