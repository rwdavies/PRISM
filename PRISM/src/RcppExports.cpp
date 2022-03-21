// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// addMatVec
Rcpp::NumericMatrix addMatVec(Rcpp::NumericMatrix A, Rcpp::NumericVector b);
RcppExport SEXP _PRISM_addMatVec(SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(addMatVec(A, b));
    return rcpp_result_gen;
END_RCPP
}
// addMatMat
Rcpp::NumericMatrix addMatMat(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B);
RcppExport SEXP _PRISM_addMatMat(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(addMatMat(A, B));
    return rcpp_result_gen;
END_RCPP
}
// increment
Rcpp::NumericVector increment(Rcpp::IntegerVector y, int yT, int xT);
RcppExport SEXP _PRISM_increment(SEXP ySEXP, SEXP yTSEXP, SEXP xTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type yT(yTSEXP);
    Rcpp::traits::input_parameter< int >::type xT(xTSEXP);
    rcpp_result_gen = Rcpp::wrap(increment(y, yT, xT));
    return rcpp_result_gen;
END_RCPP
}
// increment2N
Rcpp::NumericVector increment2N(Rcpp::NumericVector y, Rcpp::NumericVector z, int yT, int xT);
RcppExport SEXP _PRISM_increment2N(SEXP ySEXP, SEXP zSEXP, SEXP yTSEXP, SEXP xTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type yT(yTSEXP);
    Rcpp::traits::input_parameter< int >::type xT(xTSEXP);
    rcpp_result_gen = Rcpp::wrap(increment2N(y, z, yT, xT));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PRISM_addMatVec", (DL_FUNC) &_PRISM_addMatVec, 2},
    {"_PRISM_addMatMat", (DL_FUNC) &_PRISM_addMatMat, 2},
    {"_PRISM_increment", (DL_FUNC) &_PRISM_increment, 3},
    {"_PRISM_increment2N", (DL_FUNC) &_PRISM_increment2N, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_PRISM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
