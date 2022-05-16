// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fr_R_cpp
double fr_R_cpp(arma::vec x, arma::mat y, double R);
RcppExport SEXP _fourier_fr_R_cpp(SEXP xSEXP, SEXP ySEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fr_R_cpp(x, y, R));
    return rcpp_result_gen;
END_RCPP
}
// fr_R_cpp_vec
List fr_R_cpp_vec(arma::vec x, arma::mat y, double R);
RcppExport SEXP _fourier_fr_R_cpp_vec(SEXP xSEXP, SEXP ySEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fr_R_cpp_vec(x, y, R));
    return rcpp_result_gen;
END_RCPP
}
// fr_Rm_cpp
arma::vec fr_Rm_cpp(arma::vec x, arma::mat y, arma::vec R);
RcppExport SEXP _fourier_fr_Rm_cpp(SEXP xSEXP, SEXP ySEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(fr_Rm_cpp(x, y, R));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fourier_fr_R_cpp", (DL_FUNC) &_fourier_fr_R_cpp, 3},
    {"_fourier_fr_R_cpp_vec", (DL_FUNC) &_fourier_fr_R_cpp_vec, 3},
    {"_fourier_fr_Rm_cpp", (DL_FUNC) &_fourier_fr_Rm_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fourier(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}