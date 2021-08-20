// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// test_ibd
Rcpp::List test_ibd(const Eigen::MatrixXd& x, const Eigen::MatrixXd& c, const Eigen::MatrixXd& lhs, const Eigen::VectorXd& rhs, const bool approx, const int maxit, const double abstol);
RcppExport SEXP _elmulttest_test_ibd(SEXP xSEXP, SEXP cSEXP, SEXP lhsSEXP, SEXP rhsSEXP, SEXP approxSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type rhs(rhsSEXP);
    Rcpp::traits::input_parameter< const bool >::type approx(approxSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ibd(x, c, lhs, rhs, approx, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// el_pairwise
Rcpp::List el_pairwise(const Eigen::MatrixXd& x, const Eigen::MatrixXd& c, const int control, const int k, const double level, const bool interval, const std::string method, const int B, const bool approx, const int nthread, const int maxit, const double abstol);
RcppExport SEXP _elmulttest_el_pairwise(SEXP xSEXP, SEXP cSEXP, SEXP controlSEXP, SEXP kSEXP, SEXP levelSEXP, SEXP intervalSEXP, SEXP methodSEXP, SEXP BSEXP, SEXP approxSEXP, SEXP nthreadSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const int >::type control(controlSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double >::type level(levelSEXP);
    Rcpp::traits::input_parameter< const bool >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    Rcpp::traits::input_parameter< const bool >::type approx(approxSEXP);
    Rcpp::traits::input_parameter< const int >::type nthread(nthreadSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(el_pairwise(x, c, control, k, level, interval, method, B, approx, nthread, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// tt
Rcpp::List tt(const Eigen::MatrixXd& x, const Eigen::MatrixXd& c, const int control, const int k, const int maxit, const double abstol);
RcppExport SEXP _elmulttest_tt(SEXP xSEXP, SEXP cSEXP, SEXP controlSEXP, SEXP kSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const int >::type control(controlSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(tt(x, c, control, k, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}
// el_mean
Rcpp::List el_mean(const Eigen::Map<Eigen::VectorXd>& theta, const Eigen::Map<Eigen::MatrixXd>& x, const int maxit, const double abstol);
RcppExport SEXP _elmulttest_el_mean(SEXP thetaSEXP, SEXP xSEXP, SEXP maxitSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(el_mean(theta, x, maxit, abstol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_elmulttest_test_ibd", (DL_FUNC) &_elmulttest_test_ibd, 7},
    {"_elmulttest_el_pairwise", (DL_FUNC) &_elmulttest_el_pairwise, 12},
    {"_elmulttest_tt", (DL_FUNC) &_elmulttest_tt, 6},
    {"_elmulttest_el_mean", (DL_FUNC) &_elmulttest_el_mean, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_elmulttest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
