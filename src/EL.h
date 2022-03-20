#ifndef EL_H_
#define EL_H_

#include "eigen_config.h"
#include <RcppEigen.h>
#include <cmath>
#include "utils.h"

class EL2
{
private:
  const std::string type;
  Eigen::VectorXd par;
  const int maxit;
  const double tol;
  const double threshold;
  const int n;
public:
  Eigen::VectorXd l;
  double nlogLR{0};
  int iterations{1};
  bool convergence{false};

  void setup(const Eigen::Ref<const Eigen::MatrixXd>& g);

  void setup(const Eigen::Ref<const Eigen::MatrixXd>& g,
             const Eigen::Ref<const Eigen::VectorXd>& w);

  // direct evaluation
  EL2(const Eigen::Ref<const Eigen::MatrixXd>& g,
      const int maxit,
      const double tol,
      const double threshold);

  // direct evaluation (weighted)
  EL2(const Eigen::Ref<const Eigen::MatrixXd>& g,
      const Eigen::Ref<const Eigen::VectorXd>& w,
      const int maxit,
      const double tol,
      const double threshold);

  // evaluation
  EL2(const std::string method,
      const Eigen::Ref<const Eigen::VectorXd>& par0,
      const Eigen::Ref<const Eigen::MatrixXd>& x,
      const int maxit,
      const double tol,
      const double threshold);

  // evaluation (weighted)
  EL2(const std::string method,
      const Eigen::Ref<const Eigen::VectorXd>& par0,
      const Eigen::Ref<const Eigen::MatrixXd>& x,
      const Eigen::Ref<const Eigen::VectorXd>& w,
      const int maxit,
      const double tol,
      const double threshold);

  // minimization
  EL2(const std::string method,
      const Eigen::Ref<const Eigen::VectorXd>& par0,
      const Eigen::Ref<const Eigen::MatrixXd>& x,
      const Eigen::Ref<const Eigen::MatrixXd>& lhs,
      const Eigen::Ref<const Eigen::VectorXd>& rhs,
      const int maxit,
      const double tol,
      const double threshold);

  // log probability
  Eigen::ArrayXd log_prob(const Eigen::Ref<const Eigen::MatrixXd>& x,
                          const Eigen::Ref<const Eigen::ArrayXd>& w);

  // log weighted probability
  Eigen::ArrayXd log_wprob(const Eigen::Ref<const Eigen::MatrixXd>& x,
                           const Eigen::Ref<const Eigen::ArrayXd>& w);
};

class PSEUDO_LOG
{
public:
  Eigen::ArrayXd dplog;
  Eigen::ArrayXd sqrt_neg_d2plog;
  double plog_sum = 0;

  PSEUDO_LOG(Eigen::VectorXd&& x);
  PSEUDO_LOG(Eigen::VectorXd&& x, const Eigen::Ref<const Eigen::VectorXd>& w);
  static Eigen::ArrayXd plog(Eigen::VectorXd&& x);
  static double sum(Eigen::VectorXd&& x);
  static double sum(Eigen::VectorXd&& x,
                    const Eigen::Ref<const Eigen::VectorXd>& w);
};
#endif
