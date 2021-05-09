// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// el_mean
Rcpp::List el_mean(arma::rowvec theta, arma::mat x, int maxit, double abstol);
RcppExport SEXP _elmulttest_el_mean(SEXP thetaSEXP, SEXP xSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(el_mean(theta, x, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// el_mean2
Rcpp::List el_mean2(const arma::vec& theta, const arma::mat& x, const int& maxit, const double& abstol);
RcppExport SEXP _elmulttest_el_mean2(SEXP thetaSEXP, SEXP xSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double& >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(el_mean2(theta, x, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// test2sample2_cpp
Rcpp::List test2sample2_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y, double b, double alpha, int maxit, double abstol);
RcppExport SEXP _elmulttest_test2sample2_cpp(SEXP xSEXP, SEXP ySEXP, SEXP bSEXP, SEXP alphaSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(test2sample2_cpp(x, y, b, alpha, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// test2sample777_cpp
Rcpp::List test2sample777_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y, double b, double alpha, int maxit, double abstol);
RcppExport SEXP _elmulttest_test2sample777_cpp(SEXP xSEXP, SEXP ySEXP, SEXP bSEXP, SEXP alphaSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(test2sample777_cpp(x, y, b, alpha, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// test_pair
Rcpp::List test_pair(const arma::mat& x, const arma::mat& c, const std::vector<int>& pair, int maxit, double abstol);
RcppExport SEXP _elmulttest_test_pair(SEXP xSEXP, SEXP cSEXP, SEXP pairSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type pair(pairSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(test_pair(x, c, pair, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// test_ibd
Rcpp::List test_ibd(const arma::mat& x, const arma::mat& c, const arma::mat& L, const arma::vec& rhs, const int& maxit, const double& abstol);
RcppExport SEXP _elmulttest_test_ibd(SEXP xSEXP, SEXP cSEXP, SEXP LSEXP, SEXP rhsSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double& >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ibd(x, c, L, rhs, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// pairwise_ibd
Rcpp::List pairwise_ibd(const arma::mat& x, const arma::mat& c, const bool& interval, const int& B, const double& level, const int& maxit, const double& abstol);
RcppExport SEXP _elmulttest_pairwise_ibd(SEXP xSEXP, SEXP cSEXP, SEXP intervalSEXP, SEXP BSEXP, SEXP levelSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const bool& >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< const int& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const double& >::type level(levelSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double& >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwise_ibd(x, c, interval, B, level, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_elmulttest_el_mean", (DL_FUNC) &_elmulttest_el_mean, 4},
    {"_elmulttest_el_mean2", (DL_FUNC) &_elmulttest_el_mean2, 4},
    {"_elmulttest_test2sample2_cpp", (DL_FUNC) &_elmulttest_test2sample2_cpp, 6},
    {"_elmulttest_test2sample777_cpp", (DL_FUNC) &_elmulttest_test2sample777_cpp, 6},
    {"_elmulttest_test_pair", (DL_FUNC) &_elmulttest_test_pair, 5},
    {"_elmulttest_test_ibd", (DL_FUNC) &_elmulttest_test_ibd, 6},
    {"_elmulttest_pairwise_ibd", (DL_FUNC) &_elmulttest_pairwise_ibd, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_elmulttest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
