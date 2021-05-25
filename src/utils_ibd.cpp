#include "utils_ibd.h"

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat centering_ibd(arma::mat x, const arma::umat& c)
{
  // centering with nonzero elements
  x.each_col([](arma::vec& v) {
    v.elem(arma::find(v)) -= arma::mean(v.elem(arma::find(v)));
  });

  return x;
}
