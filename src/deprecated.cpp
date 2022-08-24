#include "deprecated.h"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <string>

// [[Rcpp::export]]
Rcpp::List pairwise(const Eigen::MatrixXd &x,
                    const Eigen::MatrixXd &c,
                    const int control = 0,
                    const int k = 1,
                    const double level = 0.05,
                    const bool interval = true,
                    const std::string method = "AMC",
                    const int B = 1e+04,
                    const int nthreads = 1,
                    const double th = 50,
                    const int maxit = 1e+04,
                    const double abstol = 1e-08)
{
  std::vector<std::array<int, 2>> pairs = comparison_pairs(x.cols(), control);
  const int m = pairs.size();
  std::vector<double> estimate(m);
  std::vector<double> statistic(m);
  std::vector<bool> convergence(m);
  const Eigen::VectorXd theta_hat =
      x.array().colwise().sum() / c.array().colwise().sum();
  for (int i = 0; i < m; ++i)
  {
    Rcpp::checkUserInterrupt();
    estimate[i] = theta_hat(pairs[i][0]) - theta_hat(pairs[i][1]);
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
    lhs(pairs[i][0]) = 1;
    lhs(pairs[i][1]) = -1;
    minEL pairwise_result =
        test_gbd_EL(theta_hat, x, c, lhs, Eigen::Matrix<double, 1, 1>(0),
                    th*100, maxit, abstol);
    statistic[i] = 2 * pairwise_result.nllr;
    convergence[i] = pairwise_result.convergence;
  }
  bool anyfail = std::any_of(convergence.begin(), convergence.end(),
                             [](bool v)
                             { return !v; });
  Eigen::ArrayXd bootstrap_statistics_pairwise(B);
  if (method == "AMC")
  {
    bootstrap_statistics_pairwise =
        bootstrap_statistics_pairwise_AMC(x, c, k, pairs, B, level);
  }
  else
  {
    bootstrap_statistics_pairwise = bootstrap_statistics_pairwise_NB(
        x, c, k, pairs, B, level, nthreads, th*100, maxit, abstol);
  }
  std::vector<double> adj_pvalues(m);
  for (int i = 0; i < m; ++i)
  {
    adj_pvalues[i] =
        static_cast<double>(
            (bootstrap_statistics_pairwise >= statistic[i]).count()) / B;
  }
  Rcpp::Function quantile("quantile");
  double cutoff = Rcpp::as<double>(quantile(bootstrap_statistics_pairwise,
                                            Rcpp::Named("probs") = 1 - level));
  Rcpp::List result;
  result["estimate"] = estimate;
  result["statistic"] = statistic;
  result["convergence"] = convergence;
  result["cutoff"] = cutoff;
  if (method == "NB" && anyfail)
  {
    bootstrap_statistics_pairwise =
        bootstrap_statistics_pairwise_AMC(x, c, k, pairs, B, level);
    cutoff = Rcpp::as<double>(quantile(bootstrap_statistics_pairwise,
                                       Rcpp::Named("probs") = 1 - level));
  }
  if (interval)
  {
    std::vector<double> lower(m);
    std::vector<double> upper(m);
    for (int i = 0; i < m; ++i)
    {
      Rcpp::checkUserInterrupt();
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
      lhs(pairs[i][0]) = 1;
      lhs(pairs[i][1]) = -1;
      std::array<double, 2> ci = pair_confidence_interval_gbd(
          theta_hat, x, c, lhs, th*100, estimate[i], cutoff);
      lower[i] = ci[0];
      upper[i] = ci[1];
    }
    result["lower"] = lower;
    result["upper"] = upper;
  }
  result["p.adj"] = adj_pvalues;
  result["k"] = k;
  result["level"] = level;
  result["method"] = method;
  result["B"] = bootstrap_statistics_pairwise.size();
  result.attr("class") = Rcpp::CharacterVector({"pairwise", "melt"});
  return result;
}

EL_deprecated::EL_deprecated(const Eigen::Ref<const Eigen::MatrixXd> &g,
                             const double th, const int maxit,
                             const double abstol)
{
  l = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  iterations = 0;
  convergence = false;
  while (!convergence && iterations != maxit)
  {
    PseudoLog_deprecated log_tmp(Eigen::VectorXd::Ones(g.rows()) + g * l);
    const Eigen::MatrixXd J = g.array().colwise() * log_tmp.sqrt_neg_d2plog;
    Eigen::VectorXd step =
        (J.transpose() * J)
            .ldlt()
            .solve(J.transpose() *
                   (log_tmp.dplog / log_tmp.sqrt_neg_d2plog).matrix());
    nllr = PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g.rows()) +
                                     g * (l + step));
    while (nllr < log_tmp.plog_sum)
    {
      step /= 2;
      nllr = PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g.rows()) +
                                       g * (l + step));
    }
    l += step;
    if (nllr > th)
    {
      break;
    }
    if (nllr - log_tmp.plog_sum < abstol)
    {
      convergence = true;
    }
    else
    {
      ++iterations;
    }
  }
}

PseudoLog_deprecated::PseudoLog_deprecated(Eigen::VectorXd &&x)
{
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -std::log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;

  dplog.resize(x.size());
  sqrt_neg_d2plog.resize(x.size());
  plog_sum = 0;

  for (unsigned int i = 0; i < x.size(); ++i)
  {
    if (x[i] < a0)
    {
      dplog[i] = a2 + 2 * a3 * x[i];
      sqrt_neg_d2plog[i] = a2 / 2;
      plog_sum += a1 + a2 * x[i] + a3 * x[i] * x[i];
    }
    else
    {
      dplog[i] = 1.0 / x[i];
      sqrt_neg_d2plog[i] = 1.0 / x[i];
      plog_sum += std::log(x[i]);
    }
  }
}

double PseudoLog_deprecated::sum(Eigen::VectorXd &&x)
{
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -std::log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;
  double out = 0;
  for (unsigned int i = 0; i < x.size(); ++i)
  {
    out += x[i] < a0 ? a1 + a2 * x[i] + a3 * x[i] * x[i] : std::log(x[i]);
  }
  return out;
}

Eigen::ArrayXd PseudoLog_deprecated::dp(Eigen::VectorXd &&x)
{
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = 2.0 * n;
  static const double a2 = -1.0 * n * n;
  for (unsigned int i = 0; i < x.size(); ++i)
  {
    if (x[i] < a0)
    {
      x[i] = a1 + a2 * x[i];
    }
    else
    {
      x[i] = 1.0 / x[i];
    }
  }
  return x;
}

Eigen::VectorXd linear_projection(const Eigen::Ref<const Eigen::VectorXd> &par,
                                  const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                                  const Eigen::Ref<const Eigen::VectorXd> &rhs)
{
  return par - lhs.transpose() * (lhs * lhs.transpose()).inverse() *
                   (lhs * par - rhs);
}

void linear_projection_void(Eigen::Ref<Eigen::VectorXd> par,
                            const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                            const Eigen::Ref<const Eigen::VectorXd> &rhs)
{
  par -=
      lhs.transpose() * (lhs * lhs.transpose()).inverse() * (lhs * par - rhs);
}

Eigen::MatrixXd bootstrap_sample(const Eigen::Ref<const Eigen::MatrixXd> &x,
                                 const Eigen::Ref<const Eigen::ArrayXi> &index)
{
  Eigen::MatrixXd out(x.rows(), x.cols());
  for (int i = 0; i < x.rows(); ++i)
  {
    out.row(i) = x.row(index(i));
  }
  return out;
}

Eigen::MatrixXd g_gbd(const Eigen::Ref<const Eigen::VectorXd> &par,
                      const Eigen::Ref<const Eigen::MatrixXd> &x,
                      const Eigen::Ref<const Eigen::MatrixXd> &c)
{
  return x - (c.array().rowwise() * par.array().transpose()).matrix();
}

Eigen::MatrixXd s_gbd(const Eigen::Ref<const Eigen::MatrixXd> &x,
                        const Eigen::Ref<const Eigen::MatrixXd> &c)
{
  Eigen::MatrixXd g =
      g_gbd(x.array().colwise().sum() / c.array().colwise().sum(), x, c);
  return (g.transpose() * g) / x.rows();
}

void lambda2theta_void(const Eigen::Ref<const Eigen::VectorXd> &l,
                       Eigen::Ref<Eigen::VectorXd> par,
                       const Eigen::Ref<const Eigen::MatrixXd> &g,
                       const Eigen::Ref<const Eigen::MatrixXd> &c,
                       const double gamma)
{
  par += gamma *
         ((PseudoLog_deprecated::dp(Eigen::VectorXd::Ones(g.rows()) + g * l)
               .matrix().asDiagonal() *
           c).array().colwise().sum().transpose() * l.array()).matrix();
}

Eigen::MatrixXd rmvn(const Eigen::MatrixXd &x, const int n)
{
  Eigen::MatrixXd I(n, x.cols());
  for (int j = 0; j < x.cols(); ++j)
  {
    for (int i = 0; i < n; ++i)
    {
      I(i, j) = R::rnorm(0, 1.0);
    }
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(x);
  return I * es.operatorSqrt();
}

minEL test_gbd_EL(const Eigen::Ref<const Eigen::VectorXd> &par0,
                  const Eigen::Ref<const Eigen::MatrixXd> &x,
                  const Eigen::Ref<const Eigen::MatrixXd> &c,
                  const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                  const Eigen::Ref<const Eigen::VectorXd> &rhs,
                  const double th,
                  const int maxit,
                  const double abstol)
{
  Eigen::VectorXd par = linear_projection(par0, lhs, rhs);
  Eigen::MatrixXd g = g_gbd(par, x, c);
  Eigen::VectorXd l = EL_deprecated(g, th).l;
  double f1 =
      PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g.rows()) + g * l);
  double gamma = 1.0 / (c.colwise().sum().mean());
  bool convergence = false;
  int iterations = 0;
  while (!convergence && iterations != maxit)
  {
    Eigen::VectorXd par_tmp = par;
    lambda2theta_void(l, par_tmp, g, c, gamma);
    linear_projection_void(par_tmp, lhs, rhs);
    Eigen::MatrixXd g_tmp = g_gbd(par_tmp, x, c);
    EL_deprecated eval(g_tmp, th);
    Eigen::VectorXd l_tmp = eval.l;
    if (!eval.convergence && iterations > 9)
    {
      return {par, l, f1, iterations, convergence};
    }
    double f0 = f1;
    f1 = PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g_tmp.rows()) +
                                   g_tmp * l_tmp);
    while (f0 < f1)
    {
      gamma /= 2;
      par_tmp = par;
      lambda2theta_void(l, par_tmp, g, c, gamma);
      linear_projection_void(par_tmp, lhs, rhs);
      g_tmp = g_gbd(par_tmp, x, c);
      l_tmp = EL_deprecated(g_tmp, th).l;
      if (gamma < abstol)
      {
        return {par, l, f0, iterations, convergence};
      }
      f1 = PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g_tmp.rows()) +
                                     g_tmp * l_tmp);
    }
    par = std::move(par_tmp);
    l = std::move(l_tmp);
    g = std::move(g_tmp);
    ++iterations;
    if (f0 - f1 < abstol)
    {
      convergence = true;
    }
  }
  return {par, l, f1, iterations, convergence};
}

double test_nllr(const Eigen::Ref<const Eigen::VectorXd> &par0,
                 const Eigen::Ref<const Eigen::MatrixXd> &x,
                 const Eigen::Ref<const Eigen::MatrixXd> &c,
                 const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                 const Eigen::Ref<const Eigen::VectorXd> &rhs,
                 const double th,
                 const int maxit,
                 const double abstol)
{
  Eigen::VectorXd par = linear_projection(par0, lhs, rhs);
  Eigen::MatrixXd g = g_gbd(par, x, c);
  Eigen::VectorXd l = EL_deprecated(g, th).l;
  double f1 =
      PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g.rows()) + g * l);
  double gamma = 1.0 / (c.colwise().sum().mean());
  bool convergence = false;
  int iterations = 0;
  while (!convergence && iterations != maxit)
  {
    Eigen::VectorXd par_tmp = par;
    lambda2theta_void(l, par_tmp, g, c, gamma);
    linear_projection_void(par_tmp, lhs, rhs);
    Eigen::MatrixXd g_tmp = g_gbd(par_tmp, x, c);
    EL_deprecated eval(g_tmp, th);
    Eigen::VectorXd l_tmp = eval.l;
    if (!eval.convergence && iterations > 9)
    {
      return f1;
    }
    double f0 = f1;
    f1 = PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g_tmp.rows()) +
                                   g_tmp * l_tmp);
    while (f0 < f1)
    {
      gamma /= 2;
      par_tmp = par;
      lambda2theta_void(l, par_tmp, g, c, gamma);
      linear_projection_void(par_tmp, lhs, rhs);
      g_tmp = g_gbd(par_tmp, x, c);
      l_tmp = EL_deprecated(g_tmp, th).l;
      if (gamma < abstol)
      {
        return f0;
      }
      f1 = PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g_tmp.rows()) +
                                     g_tmp * l_tmp);
    }
    par = std::move(par_tmp);
    l = std::move(l_tmp);
    g = std::move(g_tmp);
    ++iterations;
    if (f0 - f1 < abstol)
    {
      convergence = true;
    }
  }
  return f1;
}

double test_nllr(const Eigen::Ref<const Eigen::MatrixXd> &x,
                 const Eigen::Ref<const Eigen::MatrixXd> &c,
                 const Eigen::Ref<const Eigen::MatrixXd> &lhs,
                 const Eigen::Ref<const Eigen::VectorXd> &rhs,
                 const double th,
                 const int maxit,
                 const double abstol)
{
  Eigen::VectorXd par = linear_projection(
      x.array().colwise().sum() / c.array().colwise().sum(), lhs, rhs);
  Eigen::MatrixXd g = g_gbd(par, x, c);
  Eigen::VectorXd l = EL_deprecated(g, th).l;
  double f1 =
      PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g.rows()) + g * l);
  double gamma = 1.0 / (c.colwise().sum().mean());
  bool convergence = false;
  int iterations = 0;
  while (!convergence && iterations != maxit)
  {
    Eigen::VectorXd par_tmp = par;
    lambda2theta_void(l, par_tmp, g, c, gamma);
    linear_projection_void(par_tmp, lhs, rhs);
    Eigen::MatrixXd g_tmp = g_gbd(par_tmp, x, c);
    EL_deprecated eval(g_tmp, th);
    Eigen::VectorXd l_tmp = eval.l;
    if (!eval.convergence && iterations > 9)
    {
      return f1;
    }
    double f0 = f1;
    f1 = PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g_tmp.rows()) +
                                   g_tmp * l_tmp);
    while (f0 < f1)
    {
      gamma /= 2;
      par_tmp = par;
      lambda2theta_void(l, par_tmp, g, c, gamma);
      linear_projection_void(par_tmp, lhs, rhs);
      g_tmp = g_gbd(par_tmp, x, c);
      l_tmp = EL_deprecated(g_tmp, th).l;
      if (gamma < abstol)
      {
        return f0;
      }
      f1 = PseudoLog_deprecated::sum(Eigen::VectorXd::Ones(g_tmp.rows()) +
                                     g_tmp * l_tmp);
    }
    par = std::move(par_tmp);
    l = std::move(l_tmp);
    g = std::move(g_tmp);
    ++iterations;
    if (f0 - f1 < abstol)
    {
      convergence = true;
    }
  }
  return f1;
}

std::vector<std::array<int, 2>> comparison_pairs(const int p,
                                                 const int control)
{
  std::vector<std::array<int, 2>> pairs;
  if (control == 0)
  {
    pairs.reserve(p * (p - 1) / 2);
    for (int i = 0; i < p - 1; ++i)
    {
      for (int j = i + 1; j < p; ++j)
      {
        pairs.emplace_back(std::array<int, 2>{i, j});
      }
    }
  }
  else
  {
    pairs.reserve(p - 1);
    for (int i = 0; i < p; ++i)
    {
      if (i == control - 1)
      {
        continue;
      }
      pairs.emplace_back(std::array<int, 2>{i, control - 1});
    }
  }
  return pairs;
}

std::array<double, 2> pair_confidence_interval_gbd(
    const Eigen::Ref<const Eigen::VectorXd> &par0,
    const Eigen::Ref<const Eigen::MatrixXd> &x,
    const Eigen::Ref<const Eigen::MatrixXd> &c,
    const Eigen::Ref<const Eigen::MatrixXd> &lhs,
    const double th,
    const double init,
    const double cutoff)
{
  double upper_lb = init;
  double upper_size = 1;
  double upper_ub = init + upper_size;
  while (2 * test_nllr(par0, x, c, lhs, Eigen::Matrix<double, 1, 1>(upper_ub),
                       th) <=
         cutoff)
  {
    upper_lb = upper_ub;
    upper_ub += upper_size;
  }
  while (upper_ub - upper_lb > 1e-04)
  {
    if (2 * test_nllr(par0, x, c, lhs,
                      Eigen::Matrix<double, 1, 1>((upper_lb + upper_ub) / 2),
                      th) >
        cutoff)
    {
      upper_ub = (upper_lb + upper_ub) / 2;
    }
    else
    {
      upper_lb = (upper_lb + upper_ub) / 2;
    }
  }
  double lower_ub = init;
  double lower_size = 1;
  double lower_lb = init - lower_size;
  while (2 * test_nllr(par0, x, c, lhs, Eigen::Matrix<double, 1, 1>(lower_lb),
                       th) <=
         cutoff)
  {
    lower_ub = lower_lb;
    lower_lb -= lower_size;
  }
  while (lower_ub - lower_lb > 1e-04)
  {
    if (2 * test_nllr(par0, x, c, lhs,
                      Eigen::Matrix<double, 1, 1>((lower_lb + lower_ub) / 2),
                      th) >
        cutoff)
    {
      lower_lb = (lower_lb + lower_ub) / 2;
    }
    else
    {
      lower_ub = (lower_lb + lower_ub) / 2;
    }
  }
  return std::array<double, 2>{lower_ub, upper_lb};
}

Eigen::ArrayXd bootstrap_statistics_pairwise_AMC(
    const Eigen::Ref<const Eigen::MatrixXd> &x,
    const Eigen::Ref<const Eigen::MatrixXd> &c,
    const int k,
    const std::vector<std::array<int, 2>> &pairs,
    const int B,
    const double level)
{
  const Eigen::MatrixXd V_hat = s_gbd(x, c);
  const Eigen::MatrixXd U_hat = rmvn(V_hat, B);
  const int m = pairs.size();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      bootstrap_statistics(B, m);
  for (int j = 0; j < m; ++j)
  {
    Eigen::RowVectorXd R = Eigen::RowVectorXd::Zero(1, x.cols());
    R(pairs[j][0]) = 1;
    R(pairs[j][1]) = -1;
    Eigen::MatrixXd A_hat = (R.transpose() * R) / (R * V_hat * R.transpose());
    bootstrap_statistics.col(j) =
        (U_hat * A_hat * U_hat.transpose()).diagonal();
  }
  for (int b = 0; b < B; ++b)
  {
    std::sort(bootstrap_statistics.row(b).data(),
              bootstrap_statistics.row(b).data() + m);
  }
  return bootstrap_statistics.col(m - k);
}

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
    const double abstol)
{
  const int n = x.rows();
  const int p = x.cols();
  const int m = pairs.size();
  const Eigen::MatrixXd x_centered =
      x - (c.array().rowwise() *
           (x.array().colwise().sum() / c.array().colwise().sum()))
              .matrix();
  const Eigen::ArrayXXi bootstrap_index =
      Eigen::Map<Eigen::ArrayXXi, Eigen::Unaligned>(
          (Rcpp::as<std::vector<int>>(Rcpp::sample(
               Rcpp::IntegerVector(Rcpp::seq(0, n - 1)), n * B, true)))
              .data(),
          n, B);
  Eigen::ArrayXd k_bootstrap_statistics = Eigen::ArrayXd::Constant(B, NA_REAL);
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  for (int b = 0; b < B; ++b)
  {
    std::vector<double> bootstrap_statistics(m);
    for (int j = 0; j < m; ++j)
    {
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, p);
      lhs(pairs[j][0]) = 1;
      lhs(pairs[j][1]) = -1;
      bootstrap_statistics[j] =
          2 * test_nllr(bootstrap_sample(x_centered, bootstrap_index.col(b)),
                        bootstrap_sample(c, bootstrap_index.col(b)), lhs,
                        Eigen::Matrix<double, 1, 1>(0), th, maxit, abstol);
    }
    std::sort(bootstrap_statistics.begin(), bootstrap_statistics.end());
    k_bootstrap_statistics[b] = bootstrap_statistics[m - k];
  }
  return k_bootstrap_statistics;
}
