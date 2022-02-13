#ifndef EL_H_
#define EL_H_

#include "eigen_config.h"
#include <RcppEigen.h>
// #include "PSEUDO_LOG.h"
#include "utils.h"

struct minEL {
  Eigen::VectorXd par;
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;
};

class EL {
public:
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;

  EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
     const double threshold,
     const int maxit = 100,
     const double abstol = 1e-8);
};

class EL2 {
public:
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;

  EL2(const Eigen::Ref<const Eigen::VectorXd>& par,
      const Eigen::Ref<const Eigen::MatrixXd>& x,
      const std::string type,
      // const double threshold,
      const int maxit,
      const double abstol);
};

class PSEUDO_LOG {
public:
  Eigen::ArrayXd dplog;
  Eigen::ArrayXd sqrt_neg_d2plog;
  double plog_sum;

  // PSEUDO_LOG(const Eigen::Ref<const Eigen::VectorXd>& x);
  PSEUDO_LOG(Eigen::VectorXd&& x);
  static double sum(Eigen::VectorXd&& x);
  static Eigen::ArrayXd dp(Eigen::VectorXd&& x);
};

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::VectorXd>& par,
                       const Eigen::Ref<const Eigen::MatrixXd>& x);

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::VectorXd>& beta,
                     const Eigen::Ref<const Eigen::MatrixXd>& x,
                     const Eigen::Ref<const Eigen::VectorXd>& y);

Eigen::MatrixXd g_lm2(const Eigen::Ref<const Eigen::VectorXd>& par,
                      const Eigen::Ref<const Eigen::MatrixXd>& data);

Eigen::VectorXd gradient_nlogLR_lm2(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data);

minEL el_test(const Eigen::Ref<const Eigen::VectorXd>& par0,
              const Eigen::Ref<const Eigen::MatrixXd>& data,
              const std::string type,
              const Eigen::Ref<const Eigen::MatrixXd>& lhs,
              const Eigen::Ref<const Eigen::VectorXd>& rhs,
              const double threshold,
              const int maxit,
              const double abstol);
#endif
