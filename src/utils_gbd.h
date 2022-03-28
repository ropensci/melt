#ifndef EL_UTILS_gbd_H_
#define EL_UTILS_gbd_H_

#include <RcppEigen.h>

struct minEL {
  Eigen::VectorXd par;
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;
};

class EL_deprecated {
public:
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;

  EL_deprecated(const Eigen::Ref<const Eigen::MatrixXd>& g,
                const double threshold,
                const int maxit = 100,
                const double abstol = 1e-8);
};

class PSEUDO_LOG_deprecated {
public:
  Eigen::ArrayXd dplog;
  Eigen::ArrayXd sqrt_neg_d2plog;
  double plog_sum;

  PSEUDO_LOG_deprecated(Eigen::VectorXd&& x);
  static double sum(Eigen::VectorXd&& x);
  static Eigen::ArrayXd dp(Eigen::VectorXd&& x);
};

Eigen::VectorXd linear_projection(
    const Eigen::Ref<const Eigen::VectorXd>& theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs);

void linear_projection_void(
    Eigen::Ref<Eigen::VectorXd> theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs);

Eigen::MatrixXd bootstrap_sample(const Eigen::Ref<const Eigen::MatrixXd>& x,
                                 const Eigen::Ref<const Eigen::ArrayXi>& index);

Eigen::MatrixXd g_gbd(const Eigen::Ref<const Eigen::VectorXd>& theta,
                      const Eigen::Ref<const Eigen::MatrixXd>& x,
                      const Eigen::Ref<const Eigen::MatrixXd>& c);

Eigen::MatrixXd cov_gbd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                        const Eigen::Ref<const Eigen::MatrixXd>& c);

Eigen::VectorXd lambda2theta_gbd(const Eigen::Ref<const Eigen::VectorXd>& lambda,
                                 const Eigen::Ref<const Eigen::VectorXd>& theta,
                                 const Eigen::Ref<const Eigen::MatrixXd>& g,
                                 const Eigen::Ref<const Eigen::MatrixXd>& c,
                                 const double gamma);

void lambda2theta_void(
                const Eigen::Ref<const Eigen::VectorXd>& lambda,
                Eigen::Ref<Eigen::VectorXd> theta,
                const Eigen::Ref<const Eigen::MatrixXd>& g,
                const Eigen::Ref<const Eigen::MatrixXd>& c,
                const double gamma);

Eigen::VectorXd approx_lambda_gbd(
                const Eigen::Ref<const Eigen::MatrixXd>& g0,
                const Eigen::Ref<const Eigen::MatrixXd>& c,
                const Eigen::Ref<const Eigen::VectorXd>& theta0,
                const Eigen::Ref<const Eigen::VectorXd>& theta1,
                const Eigen::Ref<const Eigen::VectorXd>& lambda0);

Eigen::MatrixXd rmvn(const Eigen::MatrixXd& x, const int n);

minEL test_gbd_EL(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                  const Eigen::Ref<const Eigen::MatrixXd>& x,
                  const Eigen::Ref<const Eigen::MatrixXd>& c,
                  const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                  const Eigen::Ref<const Eigen::VectorXd>& rhs,
                  const double threshold,
                  const int maxit = 1000,
                  const double abstol = 1e-8);

double test_nlogLR(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                   const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const double threshold,
                   const int maxit = 1000,
                   const double abstol = 1e-8);

double test_nlogLR(const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const double threshold,
                   const int maxit = 1000,
                   const double abstol = 1e-8);
#endif
