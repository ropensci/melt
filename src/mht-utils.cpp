#include "mht-utils.h"

Eigen::RowVectorXd rmvn(const Eigen::Ref<const Eigen::MatrixXd>& sqrt) {
  Eigen::RowVectorXd u(sqrt.cols());
  for (int i = 0; i < sqrt.cols(); ++i) {
    u(i) = R::rnorm(0, 1.0);
  }
  return u * sqrt;
}

Eigen::MatrixXd w_mean(const Eigen::Ref<const Eigen::MatrixXd>& x)
{
  return -Eigen::MatrixXd::Identity(x.cols(), x.cols());
}

Eigen::MatrixXd dg0_inv(const std::string method,
                        const Eigen::Ref<const Eigen::MatrixXd>& x)
{
  std::map<std::string, std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd>&)>>
        w_map{{{"mean", w_mean},
               {"lm", w_mean},
               {"gaussian_identity", w_mean},
               {"gaussian_log", w_mean},
               {"gaussian_inverse", w_mean},
               {"binomial_logit", w_mean},
               {"binomial_probit", w_mean},
               {"binomial_log", w_mean},
               {"poisson_log", w_mean},
               {"poisson_identity", w_mean},
               {"poisson_sqrt", w_mean},
               {"quasibinomial_logit", w_mean}}};
  return w_map[method](x);
}





Eigen::MatrixXd cov(const std::string method,
                    const Eigen::Ref<const Eigen::VectorXd>& est,
                    const Eigen::Ref<const Eigen::MatrixXd>& x)
{
  // const Eigen::MatrixXd g = g_fn(x, est);
  // std::map<std::string, std::function<Eigen::MatrixXd(
  //     const Eigen::Ref<const Eigen::MatrixXd>&,
  //     const Eigen::Ref<const Eigen::VectorXd>&)>>
  //       g_map{{{"mean", g_mean},
  //       {"lm", g_lm},
  //       {"gaussian_identity", g_lm},
  //       {"gaussian_log", g_gauss_log},
  //       {"gaussian_inverse", g_gauss_inverse},
  //       {"binomial_logit", g_bin_logit},
  //       {"binomial_probit", g_bin_probit},
  //       {"binomial_log", g_bin_log},
  //       {"poisson_log", g_poi_log},
  //       {"poisson_identity", g_poi_identity},
  //       {"poisson_sqrt", g_poi_sqrt},
  //       {"quasibinomial_logit", g_qbin_logit}}};
  std::map<std::string, std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&)>>
        g_map{{{"mean", g_mean},
        {"lm", g_mean},
        {"gaussian_identity", g_mean},
        {"gaussian_log", g_mean},
        {"gaussian_inverse", g_mean},
        {"binomial_logit", g_mean},
        {"binomial_probit", g_mean},
        {"binomial_log", g_mean},
        {"poisson_log", g_mean},
        {"poisson_identity", g_mean},
        {"poisson_sqrt", g_mean},
        {"quasibinomial_logit", g_mean}}};
  return (1.0 / x.rows()) *
    ((g_map[method](x, est)).transpose() * g_map[method](x, est));
}

Eigen::MatrixXd ahat(const Eigen::Ref<const Eigen::MatrixXd>& j,
                     const Eigen::Ref<const Eigen::MatrixXd>& w,
                     const Eigen::Ref<const Eigen::MatrixXd>& s)
{
  const Eigen::MatrixXd tmp = j * w;
  return tmp.transpose() * ((tmp * s * tmp.transpose()).inverse()) * tmp;
}











// [[Rcpp::export]]
Eigen::MatrixXd cov2(const std::string method,
                    const Eigen::Map<Eigen::VectorXd>& est,
                    const Eigen::Map<Eigen::MatrixXd>& x)
{
  // const Eigen::MatrixXd g = g_fn(x, est);
  // std::map<std::string, std::function<Eigen::MatrixXd(
  //     const Eigen::Ref<const Eigen::MatrixXd>&,
  //     const Eigen::Ref<const Eigen::VectorXd>&)>>
  //       g_map{{{"mean", g_mean},
  //       {"lm", g_lm},
  //       {"gaussian_identity", g_lm},
  //       {"gaussian_log", g_gauss_log},
  //       {"gaussian_inverse", g_gauss_inverse},
  //       {"binomial_logit", g_bin_logit},
  //       {"binomial_probit", g_bin_probit},
  //       {"binomial_log", g_bin_log},
  //       {"poisson_log", g_poi_log},
  //       {"poisson_identity", g_poi_identity},
  //       {"poisson_sqrt", g_poi_sqrt},
  //       {"quasibinomial_logit", g_qbin_logit}}};
  std::map<std::string, std::function<Eigen::MatrixXd(
      const Eigen::Ref<const Eigen::MatrixXd>&,
      const Eigen::Ref<const Eigen::VectorXd>&)>>
        g_map{{{"mean", g_mean},
        {"lm", g_mean},
        {"gaussian_identity", g_mean},
        {"gaussian_log", g_mean},
        {"gaussian_inverse", g_mean},
        {"binomial_logit", g_mean},
        {"binomial_probit", g_mean},
        {"binomial_log", g_mean},
        {"poisson_log", g_mean},
        {"poisson_identity", g_mean},
        {"poisson_sqrt", g_mean},
        {"quasibinomial_logit", g_mean}}};
  return (1.0 / x.rows()) *
    ((g_map[method](x, est)).transpose() * g_map[method](x, est));
}
