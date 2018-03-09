#include <Rcpp.h>
#include <RcppEigen.h>
#include "matrix4.h"
#include "matrix-varia.h"
#include "wald.h"
#include <ctime>
#include <cmath>

#define scalar double
#define MATRIX MatrixXd
#define VECTOR VectorXd

using namespace Rcpp;
using namespace Eigen;

//[[Rcpp::export]]
List GxE_lm_quanti(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, 
                   NumericMatrix X, NumericMatrix I_, int beg, int end) {

  int n = Y.size();
  int r = X.ncol();

  // recopiage des matrices... en float
  MATRIX y(n,1);
  MATRIX x(n,r);
  MATRIX I(n,3);
  for(int i = 0; i < n; i++) y(i,0) = Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = X(i,j);
  
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < 3; j++)
      I(i,j) = I_(i,j);

  // Matrix used in estimation computations
  MATRIX xtx(r,r);
  MATRIX xtx_i(r,r); 
  MATRIX xty(r,1);
  MATRIX xtx_ixty(r,1);
  MATRIX xtI(r,3);
  MATRIX xtx_ixtI(r,3);
  MATRIX Ity(3,1);
  MATRIX W(3,3);
  MATRIX W_i(3,3);
  MATRIX VAR(3,3);
  MATRIX VAR2(2,2);

  scalar d, log_d;
  
  // Fixed matrix
  xtx = x.transpose() * x;
  xty = x.transpose() * y;
  sym_inverse(xtx, xtx_i, log_d, d, 1e-5);
  xtx_ixty = xtx_i * xty;

  // declare vectors containing result
  VECTOR beta(r);
  VECTOR gamma(3);
  VECTOR gamma2(2);
  scalar s2, wald;
  
  NumericVector BETA_E(end-beg+1);
  NumericVector BETA_SNP(end-beg+1);
  NumericVector BETA_ExSNP(end-beg+1);
  NumericVector SDBETA_E(end-beg+1);
  NumericVector SDBETA_SNP(end-beg+1);
  NumericVector SDBETA_ExSNP(end-beg+1);
  NumericVector COVBETA_E_SNP(end-beg+1);
  NumericVector COVBETA_E_ExSNP(end-beg+1);
  NumericVector COVBETA_SNP_ExSNP(end-beg+1);
  NumericVector WALD_3df(end-beg+1);
  NumericVector WALD_2df(end-beg+1);
  
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      BETA_E(i-beg) = NAN;
      BETA_SNP(i-beg) = NAN;
      BETA_ExSNP(i-beg) = NAN;
      SDBETA_E(i-beg) = NAN;
      SDBETA_SNP(i-beg) = NAN;
      SDBETA_ExSNP(i-beg) = NAN;
      COVBETA_E_SNP(i-beg) = NAN;
      COVBETA_E_ExSNP(i-beg) = NAN;
      COVBETA_SNP_ExSNP(i-beg) = NAN;
	  WALD_3df(i-beg) = NAN;
	  WALD_2df(i-beg) = NAN;
      continue;
    }

    // remplir SNP
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        I(4*ii+ss, 1) = ((xx&3) != 3)?(xx&3):mu(i);
        I(4*ii+ss, 2) = I(4*ii+ss, 0)*I(4*ii+ss, 1);
        xx >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        I(4*ii+ss, 1) = ((xx&3) != 3)?(xx&3):mu(i);
        I(4*ii+ss, 2) = I(4*ii+ss, 0)*I(4*ii+ss, 1);
        xx >>= 2;
      }
    }

	// Matrix depending on SNP
	Ity = I.transpose() * y;
    xtI = x.transpose() * I;
    xtx_ixtI = xtx_i * xtI;
    W = I.transpose() * I - xtI.transpose() * ( xtx_i * xtI ); 
    sym_inverse(W, W_i, log_d, d, 1e-5);
   
    beta = xtx_ixty +  xtx_ixtI * ( W_i * ( xtx_ixtI.transpose() * xty ) ) - xtx_ixtI * ( W_i * Ity );
	gamma = -W_i * ( xtx_ixtI.transpose() * xty ) + W_i * Ity;
	
	s2 = (y - x*beta - I*gamma).squaredNorm()/(n-r-3);
    
	VAR = s2 * W_i;
	
    BETA_E(i-beg) = gamma(0);
    BETA_SNP(i-beg) = gamma(1);
    BETA_ExSNP(i-beg) = gamma(2);
    SDBETA_E(i-beg) = sqrt(VAR(0,0));
    SDBETA_SNP(i-beg) = sqrt(VAR(1,1));
    SDBETA_ExSNP(i-beg) = sqrt(VAR(2,2));
    COVBETA_E_SNP(i-beg) = VAR(0,1);
    COVBETA_E_ExSNP(i-beg) = VAR(0,2);
    COVBETA_SNP_ExSNP(i-beg) = VAR(1,2);
	
    gamma2 = gamma.tail(2);
	VAR2 = VAR.bottomRightCorner(2,2);
	
	wald_compute(gamma, VAR, wald);
	WALD_3df(i-beg) = wald;
	
	wald_compute(gamma2, VAR2, wald);
	WALD_2df(i-beg) = wald;
  }
  
  List R;
  R["beta_E"] = BETA_E;
  R["beta_SNP"] = BETA_SNP;
  R["beta_ExSNP"] = BETA_ExSNP;
  R["sd_E"] = SDBETA_E;
  R["sd_SNP"] = SDBETA_SNP;
  R["sd_ExSNP"] = SDBETA_ExSNP;
  R["cov_E_SNP"] =  COVBETA_E_SNP;
  R["cov_E_ExSNP"] = COVBETA_E_ExSNP;
  R["cov_SNP_ExSNP"] = COVBETA_SNP_ExSNP;
  R["Wald_3df"] = WALD_3df;
  R["Wald_2df"] = WALD_2df;
  return R;
}


RcppExport SEXP gg_GxE_lm_quanti(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP I_SEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type I_(I_SEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    __result = Rcpp::wrap(GxE_lm_quanti(pA, mu, Y, X, I_, beg, end));
    return __result;
END_RCPP
}

