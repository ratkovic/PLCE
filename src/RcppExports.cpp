// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// basis_rcpp
Rcpp::List basis_rcpp(const arma::colvec& matvec, const arma::colvec& resvec, const arma::colvec& onesvec, const arma::colvec& theta);
RcppExport SEXP _PLCE_basis_rcpp(SEXP matvecSEXP, SEXP resvecSEXP, SEXP onesvecSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type matvec(matvecSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type resvec(resvecSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type onesvec(onesvecSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(basis_rcpp(matvec, resvec, onesvec, theta));
    return rcpp_result_gen;
END_RCPP
}
// basisvec_rcpp
Rcpp::List basisvec_rcpp(const arma::colvec& matvec, const arma::colvec& resvec, const arma::colvec& onesvec, const arma::colvec& theta);
RcppExport SEXP _PLCE_basisvec_rcpp(SEXP matvecSEXP, SEXP resvecSEXP, SEXP onesvecSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type matvec(matvecSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type resvec(resvecSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type onesvec(onesvecSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(basisvec_rcpp(matvec, resvec, onesvec, theta));
    return rcpp_result_gen;
END_RCPP
}
// corbases
Rcpp::List corbases(arma::vec treat, arma::mat X, arma::mat inter_schedule, arma::vec one_to_n);
RcppExport SEXP _PLCE_corbases(SEXP treatSEXP, SEXP XSEXP, SEXP inter_scheduleSEXP, SEXP one_to_nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type treat(treatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inter_schedule(inter_scheduleSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type one_to_n(one_to_nSEXP);
    rcpp_result_gen = Rcpp::wrap(corbases(treat, X, inter_schedule, one_to_n));
    return rcpp_result_gen;
END_RCPP
}
// bayesLasso
Rcpp::List bayesLasso(arma::vec y, arma::mat X, arma::vec alpha);
RcppExport SEXP _PLCE_bayesLasso(SEXP ySEXP, SEXP XSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(bayesLasso(y, X, alpha));
    return rcpp_result_gen;
END_RCPP
}
// fastres_cpp
List fastres_cpp(const arma::mat& X, const arma::mat& y, const arma::colvec& w);
RcppExport SEXP _PLCE_fastres_cpp(SEXP XSEXP, SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(fastres_cpp(X, y, w));
    return rcpp_result_gen;
END_RCPP
}
// rcppClamp
NumericVector rcppClamp(NumericVector x, double mi, double ma);
RcppExport SEXP _PLCE_rcppClamp(SEXP xSEXP, SEXP miSEXP, SEXP maSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mi(miSEXP);
    Rcpp::traits::input_parameter< double >::type ma(maSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppClamp(x, mi, ma));
    return rcpp_result_gen;
END_RCPP
}
// which_maxCpp
int which_maxCpp(NumericVector v);
RcppExport SEXP _PLCE_which_maxCpp(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(which_maxCpp(v));
    return rcpp_result_gen;
END_RCPP
}
// arma_var
arma::vec arma_var(const arma::vec& v1);
RcppExport SEXP _PLCE_arma_var(SEXP v1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type v1(v1SEXP);
    rcpp_result_gen = Rcpp::wrap(arma_var(v1));
    return rcpp_result_gen;
END_RCPP
}
// arma_matvar
arma::mat arma_matvar(const arma::mat& v1);
RcppExport SEXP _PLCE_arma_matvar(SEXP v1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type v1(v1SEXP);
    rcpp_result_gen = Rcpp::wrap(arma_matvar(v1));
    return rcpp_result_gen;
END_RCPP
}
// arma_ginv
arma::mat arma_ginv(const arma::mat& v1);
RcppExport SEXP _PLCE_arma_ginv(SEXP v1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type v1(v1SEXP);
    rcpp_result_gen = Rcpp::wrap(arma_ginv(v1));
    return rcpp_result_gen;
END_RCPP
}
// arma_symm
arma::mat arma_symm(const arma::mat& v1);
RcppExport SEXP _PLCE_arma_symm(SEXP v1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type v1(v1SEXP);
    rcpp_result_gen = Rcpp::wrap(arma_symm(v1));
    return rcpp_result_gen;
END_RCPP
}
// solve_cpp
Rcpp::List solve_cpp(arma::mat XpX, arma::colvec Xpy);
RcppExport SEXP _PLCE_solve_cpp(SEXP XpXSEXP, SEXP XpySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type XpX(XpXSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Xpy(XpySEXP);
    rcpp_result_gen = Rcpp::wrap(solve_cpp(XpX, Xpy));
    return rcpp_result_gen;
END_RCPP
}
// sd_cpp
double sd_cpp(arma::colvec v1);
RcppExport SEXP _PLCE_sd_cpp(SEXP v1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type v1(v1SEXP);
    rcpp_result_gen = Rcpp::wrap(sd_cpp(v1));
    return rcpp_result_gen;
END_RCPP
}
// median_cpp
double median_cpp(arma::colvec v1);
RcppExport SEXP _PLCE_median_cpp(SEXP v1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type v1(v1SEXP);
    rcpp_result_gen = Rcpp::wrap(median_cpp(v1));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PLCE_basis_rcpp", (DL_FUNC) &_PLCE_basis_rcpp, 4},
    {"_PLCE_basisvec_rcpp", (DL_FUNC) &_PLCE_basisvec_rcpp, 4},
    {"_PLCE_corbases", (DL_FUNC) &_PLCE_corbases, 4},
    {"_PLCE_bayesLasso", (DL_FUNC) &_PLCE_bayesLasso, 3},
    {"_PLCE_fastres_cpp", (DL_FUNC) &_PLCE_fastres_cpp, 3},
    {"_PLCE_rcppClamp", (DL_FUNC) &_PLCE_rcppClamp, 3},
    {"_PLCE_which_maxCpp", (DL_FUNC) &_PLCE_which_maxCpp, 1},
    {"_PLCE_arma_var", (DL_FUNC) &_PLCE_arma_var, 1},
    {"_PLCE_arma_matvar", (DL_FUNC) &_PLCE_arma_matvar, 1},
    {"_PLCE_arma_ginv", (DL_FUNC) &_PLCE_arma_ginv, 1},
    {"_PLCE_arma_symm", (DL_FUNC) &_PLCE_arma_symm, 1},
    {"_PLCE_solve_cpp", (DL_FUNC) &_PLCE_solve_cpp, 2},
    {"_PLCE_sd_cpp", (DL_FUNC) &_PLCE_sd_cpp, 1},
    {"_PLCE_median_cpp", (DL_FUNC) &_PLCE_median_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_PLCE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
