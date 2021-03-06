// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sparseADMM
Rcpp::List sparseADMM(arma::sp_mat x, arma::mat y, const Rcpp::List control);
RcppExport SEXP _solvers_sparseADMM(SEXP xSEXP, SEXP ySEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseADMM(x, y, control));
    return rcpp_result_gen;
END_RCPP
}
// denseADMM
Rcpp::List denseADMM(arma::mat x, arma::mat y, const Rcpp::List control);
RcppExport SEXP _solvers_denseADMM(SEXP xSEXP, SEXP ySEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(denseADMM(x, y, control));
    return rcpp_result_gen;
END_RCPP
}
// sparseFISTA
Rcpp::List sparseFISTA(arma::sp_mat x, arma::mat y, const Rcpp::List control);
RcppExport SEXP _solvers_sparseFISTA(SEXP xSEXP, SEXP ySEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseFISTA(x, y, control));
    return rcpp_result_gen;
END_RCPP
}
// denseFISTA
Rcpp::List denseFISTA(arma::mat x, arma::mat y, const Rcpp::List control);
RcppExport SEXP _solvers_denseFISTA(SEXP xSEXP, SEXP ySEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(denseFISTA(x, y, control));
    return rcpp_result_gen;
END_RCPP
}
// sparsePN
Rcpp::List sparsePN(arma::sp_mat x, arma::mat y, const Rcpp::List control);
RcppExport SEXP _solvers_sparsePN(SEXP xSEXP, SEXP ySEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(sparsePN(x, y, control));
    return rcpp_result_gen;
END_RCPP
}
// densePN
Rcpp::List densePN(arma::mat x, arma::mat y, const Rcpp::List control);
RcppExport SEXP _solvers_densePN(SEXP xSEXP, SEXP ySEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(densePN(x, y, control));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_solvers_sparseADMM", (DL_FUNC) &_solvers_sparseADMM, 3},
    {"_solvers_denseADMM", (DL_FUNC) &_solvers_denseADMM, 3},
    {"_solvers_sparseFISTA", (DL_FUNC) &_solvers_sparseFISTA, 3},
    {"_solvers_denseFISTA", (DL_FUNC) &_solvers_denseFISTA, 3},
    {"_solvers_sparsePN", (DL_FUNC) &_solvers_sparsePN, 3},
    {"_solvers_densePN", (DL_FUNC) &_solvers_densePN, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_solvers(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
