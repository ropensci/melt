// #ifndef MHT_UTILS_H_
// #define MHT_UTILS_H_
//
// #include <RcppEigen.h>
//
// Eigen::RowVectorXd rmvn(const Eigen::Ref<const Eigen::MatrixXd>& sqrt);
//
// Eigen::MatrixXd w_mean(const Eigen::Ref<const Eigen::MatrixXd>& x);
// Eigen::MatrixXd w_lm(const Eigen::Ref<const Eigen::MatrixXd>& x);
//
// Eigen::MatrixXd dg0_inv(const std::string method,
//                         const Eigen::Ref<const Eigen::MatrixXd>& x);
//
// Eigen::MatrixXd cov(const std::string method,
//                     const Eigen::Ref<const Eigen::VectorXd>& est,
//                     const Eigen::Ref<const Eigen::MatrixXd>& x);
//
// Eigen::MatrixXd ahat(const Eigen::Ref<const Eigen::MatrixXd>& j,
//                      const Eigen::Ref<const Eigen::MatrixXd>& w,
//                      const Eigen::Ref<const Eigen::MatrixXd>& s);
//
// double cv_mvchisq(const std::string method,
//                   const Eigen::Ref<const Eigen::VectorXd>& est,
//                   const Eigen::Ref<const Eigen::MatrixXd>& x,
//                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
//                   const Eigen::Ref<const Eigen::VectorXi>& q,
//                   const int m,
//                   const double level,
//                   const int B,
//                   const int nthreads);
// #endif
