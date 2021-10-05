#ifndef EL_UTILS_pairwise_H_
#define EL_UTILS_pairwise_H_

#include "EL.h"
#include "utils.h"
#include "utils_gbd.h"
#include "nb_progress_bar.h"
#include <progress.hpp>
#include <array>
#ifdef _OPENMP
#include <omp.h>
#endif


std::vector<std::array<int, 2>> comparison_pairs(
    const int p, const int control);

std::array<double, 2> pair_confidence_interval_gbd(
    const Eigen::Ref<const Eigen::VectorXd>& theta0,
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const double threshold,
    const double init,
    const double cutoff);

Eigen::ArrayXd bootstrap_statistics_pairwise_AMC(
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const int k,
    const std::vector<std::array<int, 2>>& pairs,
    const int B,
    const double level);

Eigen::ArrayXd bootstrap_statistics_pairwise_NB(
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const int k,
    const std::vector<std::array<int, 2>>& pairs,
    const int B,
    const double level,
    const int nthread,
    const bool progress,
    const double threshold,
    const int maxit,
    const double abstol);
#endif
