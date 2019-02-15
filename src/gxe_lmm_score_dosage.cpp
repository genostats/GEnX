#include <Rcpp.h>
#include "gxe_lmm_score.h"
#include "gaston.utils/snp_filler_dosages.h"

//[[Rcpp::export]]
List GxE_lmm_score_dosage(CharacterVector filename, NumericVector PY, NumericMatrix P, 
                  NumericVector E, int df, int beg, int end ) {
    
  snp_filler_dosages<double> S(filename, beg, end, PY.size());

  GxE_lmm_score<double> x(PY, P, E, S);
  
  if (df == 2) x.run_tests_2df();
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
  R["score"] = S.L["score"];
  
  return R;
}


RcppExport SEXP gg_GxE_lmm_score_dosage(SEXP filenameSEXP, SEXP PYSEXP, SEXP PSEXP, SEXP ESEXP, SEXP dfSEXP, SEXP begSEXP,
                                      SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_lmm_score_dosage(filename, PY, P, E, df, beg, end));
    return rcpp_result_gen;                                                                                                                                             
END_RCPP                                                                                                                                                                
}
