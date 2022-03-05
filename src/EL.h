#ifndef EL_H_
#define EL_H_

#include "eigen_config.h"
#include <RcppEigen.h>
#include <cmath>
#include "utils.h"

struct minEL {
  Eigen::VectorXd par;
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;
};

class EL
{
private:
  const int n;
public:
  Eigen::VectorXd lambda;
  double nlogLR = 0;
  int iterations = 1;
  bool convergence = false;

  EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
     const int maxit,
     const double abstol,
     const double threshold);

  EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
     const Eigen::Ref<const Eigen::ArrayXd>& w,
     const int maxit,
     const double abstol,
     const double threshold);

  // // log probability
  // Eigen::ArrayXd log_prob(const Eigen::Ref<const Eigen::MatrixXd>& g,
  //                         const Eigen::Ref<const Eigen::ArrayXd>& w) const;
};

class EL2
{
private:
  const std::string type;
  Eigen::VectorXd par;
public:
  Eigen::VectorXd lambda;
  double nlogLR = 0;
  int iterations = 1;
  bool convergence = false;

  // evaluation
  EL2(const Eigen::Ref<const Eigen::VectorXd>& par0,
      const Eigen::Ref<const Eigen::MatrixXd>& x,
      const std::string method,
      const int maxit,
      const double abstol,
      const double threshold);

  // evaluation (weighted)
  EL2(const Eigen::Ref<const Eigen::VectorXd>& par0,
      const Eigen::Ref<const Eigen::MatrixXd>& x,
      const Eigen::Ref<const Eigen::ArrayXd>& w,
      const std::string method,
      const int maxit,
      const double abstol,
      const double threshold);

  // minimization
  EL2(const Eigen::Ref<const Eigen::VectorXd>& par0,
      const Eigen::Ref<const Eigen::MatrixXd>& x,
      const std::string method,
      const Eigen::Ref<const Eigen::MatrixXd>& lhs,
      const Eigen::Ref<const Eigen::VectorXd>& rhs,
      const int maxit,
      const double abstol,
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

  // PSEUDO_LOG(const Eigen::Ref<const Eigen::VectorXd>& x);
  PSEUDO_LOG(Eigen::VectorXd&& x);
  PSEUDO_LOG(Eigen::VectorXd&& x, Eigen::ArrayXd&& w);
  static Eigen::ArrayXd plog(Eigen::ArrayXd&& x);
  static double sum(Eigen::VectorXd&& x);
  static double sum(Eigen::VectorXd&& x, Eigen::ArrayXd&& w);
  static Eigen::ArrayXd dp(Eigen::VectorXd&& x);
};

double th_nlogLR(const int p, const Rcpp::Nullable<double> threshold);

Eigen::MatrixXd g_mean(const Eigen::Ref<const Eigen::VectorXd>& par,
                       const Eigen::Ref<const Eigen::MatrixXd>& x);

Eigen::MatrixXd g_lm(const Eigen::Ref<const Eigen::VectorXd>& par,
                     const Eigen::Ref<const Eigen::MatrixXd>& data);

Eigen::VectorXd gr_nlogLR_lm(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& data);
#endif
