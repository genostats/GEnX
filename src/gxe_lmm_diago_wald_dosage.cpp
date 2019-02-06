#include <Rcpp.h>
#include "gxe_lmm_diago_wald.h"
#include "gaston.utils/snp_filler_dosages.h"


//[[Rcpp::export]]
List GxE_lmm_wald_dosage(CharacterVector filename, NumericVector Y, NumericMatrix X, 
                  int p, NumericVector Sigma, NumericMatrix U, int df, int beg, int end, double tol ) {
    
  snp_filler_dosages<double> S(filename, beg, end, Y.size());

  GxE_lmm_wald x(Y, X, p, Sigma, U, tol, S);
  
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
  R["beta_E"] = S.L["beta_E"];
  R["beta_SNP"] = S.L["beta_SNP"];
  R["beta_ExSNP"] = S.L["beta_ExSNP"];
  R["sd_E"] = S.L["sd_E"];
  R["sd_SNP"] = S.L["sd_SNP"];
  R["sd_ExSNP"] = S.L["sd_ExSNP"];
  R["cov_E_SNP"] =  S.L["cov_E_SNP"];
  R["cov_E_ExSNP"] = S.L["cov_E_ExSNP"];
  R["cov_SNP_ExSNP"] = S.L["cov_SNP_ExSNP"];
  R["Wald"] = S.L["Wald"];
  
  return R;
}


RcppExport SEXP gg_GxE_lmm_wald_dosage(SEXP filenameSEXP, SEXP YSEXP, SEXP XSEXP, SEXP pSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP dfSEXP, SEXP begSEXP,
                                      SEXP endSEXP, SEXP tolSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(GxE_lmm_wald_dosage(filename, Y, X, p, Sigma, U, df, beg, end, tol));
    return rcpp_result_gen;                                                                                                                                             
END_RCPP                                                                                                                                                                
}
