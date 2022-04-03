#ifndef EL_H_
#define EL_H_

#include "eigen_config.h"
#include <RcppEigen.h>
#include "utils.h"

class EL
{
public:
  // members
  Eigen::VectorXd l;  // Lagrange multiplier
  double nllr{0};     // negative log likelihood ratio
  int iter{0};        // iterations performed in optimization
  bool conv{false};   // convergence status

  // constructors
  // direct evaluation
  EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
     const int maxit,
     const double tol,
     const double th);
  // direct evaluation (weighted)
  EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
     const Eigen::Ref<const Eigen::ArrayXd>& w,
     const int maxit,
     const double tol,
     const double th);
  // evaluation
  EL(const std::string method,
     const Eigen::Ref<const Eigen::VectorXd>& par0,
     const Eigen::Ref<const Eigen::MatrixXd>& x,
     const int maxit,
     const double tol,
     const double th);
  // evaluation (weighted)
  EL(const std::string method,
     const Eigen::Ref<const Eigen::VectorXd>& par0,
     const Eigen::Ref<const Eigen::MatrixXd>& x,
     const Eigen::Ref<const Eigen::ArrayXd>& w,
     const int maxit,
     const double tol,
     const double th);

  // functions for constructors
  void set_el(const Eigen::Ref<const Eigen::MatrixXd>& g);
  void set_el(const Eigen::Ref<const Eigen::MatrixXd>& g,
              const Eigen::Ref<const Eigen::ArrayXd>& w);
  std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd>&,
                                const Eigen::Ref<const Eigen::VectorXd>&)>
    set_g_fcn(const std::string method);

  // methods
  // log probability
  Eigen::ArrayXd logp(const Eigen::Ref<const Eigen::MatrixXd>& x) const;
  Eigen::ArrayXd logp(const Eigen::Ref<const Eigen::MatrixXd>& x,
                      const Eigen::Ref<const Eigen::ArrayXd>& w) const;
  Eigen::ArrayXd logp_g(const Eigen::Ref<const Eigen::MatrixXd>& g) const;
  Eigen::ArrayXd logp_g(const Eigen::Ref<const Eigen::MatrixXd>& g,
                        const Eigen::Ref<const Eigen::ArrayXd>& w) const;
  double loglik() const;
  double loglik(const Eigen::Ref<const Eigen::ArrayXd>& w) const;

private:
  // members
  const Eigen::VectorXd par;  // parameter value
  const int maxit;            // maximum number of iterations
  const double tol;           // relative convergence tolerance
  const double th;            // threshold value for -logLR
  const int n;                // sample size
  // estimating function
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&)> g_fcn;
};

class MINEL
{
public:
  // members
  Eigen::VectorXd par;  // parameter value (solution of optimization problem)
  Eigen::VectorXd l;    // Lagrange multiplier
  double nllr{0};       // negative log likelihood ratio
  int iter{0};          // iterations performed in optimization
  bool conv{false};     // convergence status

  // constructors
  // minimization
  MINEL(const std::string method,
        const Eigen::Ref<const Eigen::VectorXd>& par0,
        const Eigen::Ref<const Eigen::MatrixXd>& x,
        const Eigen::Ref<const Eigen::MatrixXd>& lhs,
        const Eigen::Ref<const Eigen::VectorXd>& rhs,
        const int maxit,
        const double tol,
        const double th);
  // minimization (weighted)
  MINEL(const std::string method,
        const Eigen::Ref<const Eigen::VectorXd>& par0,
        const Eigen::Ref<const Eigen::MatrixXd>& x,
        const Eigen::Ref<const Eigen::ArrayXd>& w,
        const Eigen::Ref<const Eigen::MatrixXd>& lhs,
        const Eigen::Ref<const Eigen::VectorXd>& rhs,
        const int maxit,
        const double tol,
        const double th);

  // functions for constructors
  std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd>&,
                                const Eigen::Ref<const Eigen::VectorXd>&)>
    set_g_fcn(const std::string method);
  std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::VectorXd>&,
                                const Eigen::Ref<const Eigen::MatrixXd>&,
                                const Eigen::Ref<const Eigen::MatrixXd>&,
                                const Eigen::Ref<const Eigen::VectorXd>&)>
    set_gr_fcn(const std::string method);
  std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::VectorXd>&,
                                const Eigen::Ref<const Eigen::MatrixXd>&,
                                const Eigen::Ref<const Eigen::MatrixXd>&,
                                const Eigen::Ref<const Eigen::ArrayXd>&)>
    set_wgr_fcn(const std::string method);

  // methods
  // log probability
  Eigen::ArrayXd logp(const Eigen::Ref<const Eigen::MatrixXd>& x) const;
  Eigen::ArrayXd logp(const Eigen::Ref<const Eigen::MatrixXd>& x,
                      const Eigen::Ref<const Eigen::ArrayXd>& w) const;
  double loglik() const;
  double loglik(const Eigen::Ref<const Eigen::ArrayXd>& w) const;

private:
  // members
  const int maxit;  // maximum number of iterations
  const double tol; // relative convergence tolerance
  const double th;  // threshold value for -logLR
  const int n;      // sample size
  // estimating function
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&)> g_fcn;
  // gradient function of negative log likelihood ratio function
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::VectorXd>&,
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&)> gr_fcn;
  const std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::VectorXd>&,
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::ArrayXd>&)> wgr_fcn;
};


class PSEUDO_LOG
{
public:
  // members
  Eigen::ArrayXd dplog;
  Eigen::ArrayXd sqrt_neg_d2plog;
  double plog_sum{0};

  // constructors
  PSEUDO_LOG(Eigen::VectorXd&& x);
  PSEUDO_LOG(Eigen::VectorXd&& x, const Eigen::Ref<const Eigen::ArrayXd>& w);

  // methods
  static Eigen::ArrayXd plog(Eigen::VectorXd&& x);
  static Eigen::ArrayXd plog(Eigen::VectorXd&& x,
                             const Eigen::Ref<const Eigen::ArrayXd>& w);
  static double sum(Eigen::VectorXd&& x);
  static double sum(Eigen::VectorXd&& x,
                    const Eigen::Ref<const Eigen::ArrayXd>& w);
};
#endif
