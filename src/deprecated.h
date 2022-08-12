#ifndef EL_deprecated_H_
#define EL_deprecated_H_

#include <RcppEigen.h>
#include <array>
#include <vector>

struct minEL
{
        Eigen::VectorXd par;
        Eigen::VectorXd l;
        double nllr;
        int iterations;
        bool convergence;
};

class EL_deprecated
{
public:
        Eigen::VectorXd l;
        double nllr;
        int iterations;
        bool convergence;

        EL_deprecated(const Eigen::Ref<const Eigen::MatrixXd> &g,
                      const double th, const int maxit = 100,
                      const double abstol = 1e-08);
};

class PseudoLog_deprecated
{
public:
        Eigen::ArrayXd dplog;
        Eigen::ArrayXd sqrt_neg_d2plog;
        double plog_sum;

        PseudoLog_deprecated(Eigen::VectorXd &&x);
        static double sum(Eigen::VectorXd &&x);
        static Eigen::ArrayXd dp(Eigen::VectorXd &&x);
};

Eigen::VectorXd linear_projection(const Eigen::Ref<const Eigen::VectorXd> &par,
                                  const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                                  const Eigen::Ref<const Eigen::VectorXd> &rhs);

void linear_projection_void(Eigen::Ref<Eigen::VectorXd> par,
                            const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                            const Eigen::Ref<const Eigen::VectorXd> &rhs);

Eigen::MatrixXd bootstrap_sample(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                 const Eigen::Ref<const Eigen::ArrayXi> &index);

Eigen::MatrixXd g_gbd(const Eigen::Ref<const Eigen::VectorXd> &par,
                      const Eigen::Ref<const Eigen::MatrixXd> &x,
                      const Eigen::Ref<const Eigen::MatrixXd> &c);

Eigen::MatrixXd s_gbd(const Eigen::Ref<const Eigen::MatrixXd> &x,
                      const Eigen::Ref<const Eigen::MatrixXd> &c);

Eigen::MatrixXd rmvn(const Eigen::MatrixXd &x, const int n);

minEL test_gbd_EL(const Eigen::Ref<const Eigen::VectorXd> &par0,
                  const Eigen::Ref<const Eigen::MatrixXd> &x,
                  const Eigen::Ref<const Eigen::MatrixXd> &c,
                  const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                  const Eigen::Ref<const Eigen::VectorXd> &rhs,
                  const double th,
                  const int maxit = 1000,
                  const double abstol = 1e-08);

double test_nllr(const Eigen::Ref<const Eigen::VectorXd> &par0,
                 const Eigen::Ref<const Eigen::MatrixXd> &x,
                 const Eigen::Ref<const Eigen::MatrixXd> &c,
                 const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                 const Eigen::Ref<const Eigen::VectorXd> &rhs,
                 const double th,
                 const int maxit = 1000,
                 const double abstol = 1e-08);

double test_nllr(const Eigen::Ref<const Eigen::MatrixXd> &x,
                 const Eigen::Ref<const Eigen::MatrixXd> &c,
                 const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                 const Eigen::Ref<const Eigen::VectorXd> &rhs,
                 const double th,
                 const int maxit = 1000,
                 const double abstol = 1e-08);

std::vector<std::array<int, 2>> comparison_pairs(const int p,
                                                 const int control);

std::array<double, 2> pair_confidence_interval_gbd(
                const Eigen::Ref<const Eigen::VectorXd> &par0,
                const Eigen::Ref<const Eigen::MatrixXd> &x,
                const Eigen::Ref<const Eigen::MatrixXd> &c,
                const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                const double th,
                const double init,
                const double cutoff);

Eigen::ArrayXd bootstrap_statistics_pairwise_AMC(
                const Eigen::Ref<const Eigen::MatrixXd> &x,
                const Eigen::Ref<const Eigen::MatrixXd> &c,
                const int k,
                const std::vector<std::array<int, 2>> &pairs,
                const int B,
                const double level);

Eigen::ArrayXd bootstrap_statistics_pairwise_NB(
    const Eigen::Ref<const Eigen::MatrixXd> &x,
    const Eigen::Ref<const Eigen::MatrixXd> &c,
    const int k,
    const std::vector<std::array<int, 2>> &pairs,
    const int B,
    const double level,
    const int nthreads,
    const double th,
    const int maxit,
    const double abstol);
#endif
