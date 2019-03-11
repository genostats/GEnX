#include <Rcpp.h>
#include "gxe_logit_wald.h"
#include "gaston.utils/snp_filler_bed.h"


//[[Rcpp::export]]
List GxE_logit_wald_bed(XPtr<matrix4> pA, NumericVector pp, NumericVector Y, NumericMatrix X, 
                  int df, int beg, int end, double tol) {
  
  snp_filler_additive_bed<double> S(pA, pp, beg, end);
  
  GxE_logit_wald x(Y, X, tol, S);
  
  if(df == 1) x.run_tests_1df();
  else if (df == 2) x.run_tests_2df();
  else if (df == 3) x.run_tests_3df();
  else stop("Bad df");
  
  return S.L;
}


RcppExport SEXP gg_GxE_logit_wald_bed(SEXP pASEXP, SEXP ppSEXP, SEXP YSEXP, SEXP XSEXP, SEXP dfSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP );
    Rcpp::traits::input_parameter< NumericVector >::type pp(ppSEXP );
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_logit_wald_bed(pA, pp, Y, X, df, beg, end, tol));
    return rcpp_result_gen;                                                                                                                                             
END_RCPP                                                                                                                                                                
}
