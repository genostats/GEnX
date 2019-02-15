#include <Rcpp.h>
#include "gxe_lmm_score.h"
#include "gaston.utils/snp_filler_bed.h"

//[[Rcpp::export]]
List GxE_lmm_score_bed(XPtr<matrix4> pA, NumericVector pp, NumericVector PY, NumericMatrix P, 
                  NumericVector E, int df, int beg, int end ) {
  
  snp_filler_additive_bed<double> S(pA, pp, beg, end);
  
  GxE_lmm_score<double> x(PY, P, E, S);
  
  if (df == 2) x.run_tests_2df();
  else if (df == 3) x.run_tests_3df();
  else stop("Bad df");
  
  return S.L;
}


RcppExport SEXP gg_GxE_lmm_score_bed(SEXP pASEXP, SEXP ppSEXP, SEXP PYSEXP, SEXP PSEXP, SEXP ESEXP, SEXP dfSEXP,
                                      SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP );
    Rcpp::traits::input_parameter< NumericVector >::type pp(ppSEXP );
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_lmm_score_bed(pA, pp, PY, P, E, df, beg, end));
    return rcpp_result_gen;                                                                                                                                             
END_RCPP                                                                                                                                                                
}
