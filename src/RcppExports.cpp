// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_prevalence_cpp
Rcpp::List get_prevalence_cpp(Rcpp::List args_params);
RcppExport SEXP _DRpower_get_prevalence_cpp(SEXP args_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args_params(args_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_prevalence_cpp(args_params));
    return rcpp_result_gen;
END_RCPP
}
// get_ICC_cpp
Rcpp::List get_ICC_cpp(Rcpp::List args_params);
RcppExport SEXP _DRpower_get_ICC_cpp(SEXP args_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args_params(args_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ICC_cpp(args_params));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DRpower_get_prevalence_cpp", (DL_FUNC) &_DRpower_get_prevalence_cpp, 1},
    {"_DRpower_get_ICC_cpp", (DL_FUNC) &_DRpower_get_ICC_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_DRpower(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
