#include "PSEUDO_LOG.h"

PSEUDO_LOG::PSEUDO_LOG(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -std::log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;

  dplog.resize(x.size());
  sqrt_neg_d2plog.resize(x.size());
  plog_sum = 0;

  for (unsigned int i = 0; i < x.size(); ++i) {
    if (x[i] < a0) {
      dplog[i] = a2 + 2 * a3 * x[i];
      sqrt_neg_d2plog[i] = a2 / 2;
      plog_sum += a1 + a2 * x[i] + a3 * x[i] * x[i];
    } else {
      dplog[i] = 1.0 / x[i];
      sqrt_neg_d2plog[i] = 1.0 / x[i];
      plog_sum += std::log(x[i]);
    }
  }
}

double PSEUDO_LOG::sum(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -std::log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;
  double out = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    out += x[i] < a0 ? a1 + a2 * x[i] + a3 * x[i] * x[i] : std::log(x[i]);
  }
  return out;
}

Eigen::ArrayXd PSEUDO_LOG::dp(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = 2.0 * n;
  static const double a2 = -1.0 * n * n;
  // Eigen::ArrayXd out(n);
  for (unsigned int i = 0; i < x.size(); ++i) {
    if (x[i] < a0) {
      x[i] = a1 + a2 * x[i];
    } else {
      x[i] = 1.0 / x[i];
    }
  }
  return x;
}


