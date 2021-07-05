#ifndef EL_UTILS_H_
#define EL_UTILS_H_

#include "eigen_config.h"
#include <RcppEigen.h>


class TEST {
public:
  double nlogLR;
  Eigen::VectorXd lambda;
  int iterations;
  bool convergence;
  TEST (const Eigen::Ref<const Eigen::MatrixXd>& g,
        const int maxit = 100,
        const double abstol = 1e-8);
};
struct EL {
  double nlogLR;
  Eigen::VectorXd lambda;
  int iterations;
  bool convergence;
};
struct minEL {
  Eigen::VectorXd theta;
  Eigen::VectorXd lambda;
  double nlogLR;
  // double p_value;
  int iterations;
  bool convergence;
};
// struct PSEUDO_LOG {
//   double plog_sum;
//   Eigen::ArrayXd dplog;
//   Eigen::ArrayXd sqrt_neg_d2plog;
// };


class PSEUDO_LOG {
public:
  double plog_sum;
  Eigen::ArrayXd dplog;
  Eigen::ArrayXd sqrt_neg_d2plog;
  // Constructor with parameters
  PSEUDO_LOG (const Eigen::Ref<const Eigen::VectorXd>& x);
};
double plog_sum(Eigen::VectorXd&& x);
Eigen::ArrayXd dplog(Eigen::VectorXd&& x);


EL getEL(const Eigen::Ref<const Eigen::MatrixXd>& g,
         const int maxit = 100,
         const double abstol = 1e-8);

std::vector<std::array<int, 2>> all_pairs(const int p);
// Eigen::VectorXd linear_projection(const Eigen::Ref<const Eigen::VectorXd>& theta,
//                                   const Eigen::Ref<const Eigen::MatrixXd>& L,
//                                   const Eigen::Ref<const Eigen::VectorXd>& rhs);


Eigen::VectorXd linear_projection(
    const Eigen::Ref<const Eigen::VectorXd>& theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs);
Eigen::VectorXd linear_projection_rref(
    Eigen::VectorXd&& theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs);

// void linear_projection2(Eigen::Ref<Eigen::VectorXd> theta,
//                         const Eigen::Ref<const Eigen::MatrixXd>& lhs,
//                         const Eigen::Ref<const Eigen::VectorXd>& rhs);

Eigen::MatrixXd bootstrap_sample(const Eigen::Ref<const Eigen::MatrixXd>& x,
                                 const Eigen::Ref<const Eigen::ArrayXi>& index);





Eigen::RowVectorXd rowset(const Eigen::Ref<const Eigen::RowVectorXd>& x,
                          const Eigen::Ref<const Eigen::RowVectorXd>& c);
Eigen::MatrixXd rowsetmat(const Eigen::Ref<const Eigen::MatrixXd>& x,
                          const Eigen::Ref<const Eigen::MatrixXd>& c);
#endif
