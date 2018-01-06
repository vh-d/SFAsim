// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rtnorm_cpp
NumericVector rtnorm_cpp(int n, double mean, double sd, double lower, double upper);
RcppExport SEXP _SFAsim_rtnorm_cpp(SEXP nSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(rtnorm_cpp(n, mean, sd, lower, upper));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SFAsim_rtnorm_cpp", (DL_FUNC) &_SFAsim_rtnorm_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_SFAsim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
