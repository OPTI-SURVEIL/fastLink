// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// calcPWDcpp
NumericMatrix calcPWDcpp(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _fastLink_calcPWDcpp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(calcPWDcpp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// m_func_par
std::vector< std::vector<arma::vec> > m_func_par(const std::vector< std::vector< std::vector<arma::vec> > > temp, 
                                                 const std::vector< std::vector< std::vector<arma::vec> > > ptemp, 
                                                 const std::vector< std::vector<arma::vec> > natemp, 
                                                 const arma::vec limit1,const arma::vec limit2, 
                                                 const arma::vec nlim1, const arma::vec nlim2, 
                                                 const arma::mat ind, const arma::vec listid, 
                                                 const std::vector< std::vector<bool>> identical,
                                                 const bool dedupe, const bool matchesLink, 
                                                 const int threads);
RcppExport SEXP _fastLink_m_func_par(SEXP tempSEXP, SEXP ptempSEXP, SEXP natempSEXP, 
                                     SEXP limit1SEXP, SEXP limit2SEXP, SEXP nlim1SEXP, 
                                     SEXP nlim2SEXP, SEXP indSEXP, SEXP listidSEXP, 
                                     SEXP identicalSEXP, SEXP dedupeSEXP,
                                     SEXP matchesLinkSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< std::vector< std::vector<arma::vec> > > >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< std::vector<arma::vec> > > >::type ptemp(ptempSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector<arma::vec> > >::type natemp(natempSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type limit1(limit1SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type limit2(limit2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type nlim1(nlim1SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type nlim2(nlim2SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type ind(indSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type listid(listidSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector<bool> > >::type identical(identicalSEXP);
    Rcpp::traits::input_parameter< const bool >::type dedupe(dedupeSEXP);
    Rcpp::traits::input_parameter< const bool >::type matchesLink(matchesLinkSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(m_func_par(temp, ptemp, natemp, limit1, limit2, nlim1, 
                                            nlim2, ind, listid, identical, dedupe,
                                            matchesLink, threads));
    return rcpp_result_gen;
END_RCPP
}
