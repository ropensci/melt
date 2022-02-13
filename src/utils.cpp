#include "utils.h"

// std::vector<std::array<int, 2>> comparison_pairs(
//     const int p, const int control) {
//   // initialize a vector of vectors
//   std::vector<std::array<int, 2>> pairs;
//   if (control == 0){
//     // the size of vector is p choose 2
//     pairs.reserve(p * (p - 1) / 2);
//     // fill in each elements(pairs)
//     for (int i = 0; i < p - 1; ++i) {
//       for (int j = i + 1; j < p; ++j) {
//         pairs.emplace_back(std::array<int, 2>{i, j});
//
//       }
//     }
//   } else {
//     // the size of vector is p - 1
//     pairs.reserve(p - 1);
//     // fill in each elements(pairs)
//     for (int i = 0; i < p; ++i) {
//       if (i == control - 1) {
//         continue;
//       }
//       pairs.emplace_back(std::array<int, 2>{i, control - 1});
//     }
//   }
//   return pairs;
// }

Eigen::VectorXd linear_projection(
    const Eigen::Ref<const Eigen::VectorXd>& theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs) {
  return theta -
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * (lhs * theta - rhs);
}
void linear_projection_void(
    Eigen::Ref<Eigen::VectorXd> theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs) {
  theta -=
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * (lhs * theta - rhs);
}

Eigen::MatrixXd bootstrap_sample(
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::ArrayXi>& index) {
  Eigen::MatrixXd out(x.rows(), x.cols());
  for (int i = 0; i < x.rows(); ++i) {
    out.row(i) = x.row(index(i));
  }
  return out;
}
