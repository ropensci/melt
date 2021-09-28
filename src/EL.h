#ifndef EL_H_
#define EL_H_

#include "eigen_config.h"
#include <RcppEigen.h>
#include "PSEUDO_LOG.h"

struct EL {
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;
};

struct minEL {
  Eigen::VectorXd theta;
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;
};

class EL2 {
public:
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;

  EL2(const Eigen::Ref<const Eigen::MatrixXd>& g,
      const double threshold,
      const int maxit = 100,
      const double abstol = 1e-8);
};

EL getEL(const Eigen::Ref<const Eigen::MatrixXd>& g,
         const int maxit = 100,
         const double abstol = 1e-8);
#endif
