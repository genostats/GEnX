#include <Rcpp.h>
#include "gxe_lmm_diago_lrt.h"
#include "gaston.utils/snp_filler_dosages.h"


//[[Rcpp::export]]
List GxE_lmm_lrt_dosage(CharacterVector filename, NumericVector Y, NumericMatrix X, 
                  int p, NumericVector Sigma, NumericMatrix U, int df, int beg, int end, double tol ) {
    
  snp_filler_dosages<double> S(filename, beg, end, Y.size());

  GxE_lmm_lrt x(Y, X, p, Sigma, U, tol, S);
  
  if(df == 1) x.run_tests_1df();
  else if (df == 2) x.run_tests_2df();
  else if (df == 3) x.run_tests_3df();
  else stop("Bad df");
  
  List R;
  
  R["id"] = wrap(S.SNP_ID);
  R["chr"] = wrap(S.CHR);
  R["pos"] = wrap(S.POS);
  R["A1"] = wrap(S.AL1);
  R["A2"] = wrap(S.AL2);
  R["freq1"] = wrap(S.F1);
  R["freq2"] = wrap(S.F2);
  R["h2"] = S.L["h2"];
  R["LRT"] = S.L["LRT"];
  
  return R;
}


RcppExport SEXP gg_GxE_lmm_lrt_dosage(SEXP filenameSEXP, SEXP YSEXP, SEXP XSEXP, SEXP pSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP dfSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);                                                                                                      
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);  
    Rcpp::traits::input_parameter< int >::type p(pSEXP);                                                                                                                
    Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP);                                                                                              
    Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP); 
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);                                                                                                              
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);                                                                                                            
    Rcpp::traits::input_parameter< int >::type end(endSEXP);                                                                                                            
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_lmm_lrt_dosage(filename, Y, X, p, Sigma, U, df, beg, end, tol));                                                                   
    return rcpp_result_gen;                                                                                                                                             
END_RCPP                                                                                                                                                                
}
