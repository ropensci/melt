#ifndef EL_H_
#define EL_H_

#include "eigen_config.h"
#include <RcppEigen.h>
#include "PSEUDO_LOG.h"
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
#endif
