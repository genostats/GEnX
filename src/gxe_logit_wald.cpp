#include <Rcpp.h>
#include "gaston/logit_model.h"
#include "gaston/matrix4.h"
#include "wald.h"
#include <ctime>
#include <cmath>
#include <iostream>
#define BLOCK 20 


// Ne marche pas avec les float au niveau de l'inversion de XtWX dans logit.h
// Je ne sais pas pourquoi...

//[[Rcpp::export]]
List GxE_logit_wald_1df(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                       int beg, int end, double tol) {
  int n = Y.size();
  int r = X.ncol();

  // recopiage des matrices... en float
  /*
  MatrixXf y(n,1);
  MatrixXf x(n,r);
  for(int i = 0; i < n; i++) y(i,0) = (float) Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = (float) X(i,j);
  */
  
  // recopiage des matrices
  MatrixXd y(n,1);
  MatrixXd x(n,r);
  for(int i = 0; i < n; i++) y(i,0) = Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = X(i,j);


  // declare vectors containing result
  NumericMatrix BETA(end-beg+1,r-3);
  NumericVector BETA_E(end-beg+1);
  NumericVector BETA_SNP(end-beg+1);
  NumericVector BETA_ExSNP(end-beg+1);
  NumericVector SDBETA_E(end-beg+1);
  NumericVector SDBETA_SNP(end-beg+1);
  NumericVector SDBETA_ExSNP(end-beg+1);
  NumericVector COVBETA_E_SNP(end-beg+1);
  NumericVector COVBETA_E_ExSNP(end-beg+1);
  NumericVector COVBETA_SNP_ExSNP(end-beg+1);
  
  //VectorXf Beta(r);
  //Beta.setZero();
  VectorXd beta(r);
  beta.setZero();
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
      continue;
    }
    // remplir les dernières colonnes de x par génotype au SNP et l'interaction(manquant -> mu)
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        x(4*ii+ss, r-2) = ((xx&3) != 3)?(xx&3):mu(i);
		x(4*ii+ss, r-1) = x(4*ii+ss, r-3)*x(4*ii+ss, r-2);
        xx >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        x(4*ii+ss, r-2) = ((xx&3) != 3)?(xx&3):mu(i);
		x(4*ii+ss, r-1) = x(4*ii+ss, r-3)*x(4*ii+ss, r-2);
        xx >>= 2;
      }
    }

	MatrixXd varbeta(r,r);
	logistic_model2<double>(y, x, beta, varbeta, tol);
	
	for(int k = 0; k < r-3; k++) {
	  BETA(i-beg,k) = beta(k);
	}
	
    BETA_E(i-beg) = beta(r-3);
    BETA_SNP(i-beg) = beta(r-2);
    BETA_ExSNP(i-beg) = beta(r-1);
    SDBETA_E(i-beg) = sqrt(varbeta(r-3,r-3));
    SDBETA_SNP(i-beg) = sqrt(varbeta(r-2,r-2));
    SDBETA_ExSNP(i-beg) = sqrt(varbeta(r-1,r-1));
    COVBETA_E_SNP(i-beg) = varbeta(r-3,r-2);
    COVBETA_E_ExSNP(i-beg) = varbeta(r-3,r-1);
    COVBETA_SNP_ExSNP(i-beg) = varbeta(r-2,r-1);
  }

  List R;
  R["beta"] = BETA;
  R["beta_E"] = BETA_E;
  R["beta_SNP"] = BETA_SNP;
  R["beta_ExSNP"] = BETA_ExSNP;
  R["sd_E"] = SDBETA_E;
  R["sd_SNP"] = SDBETA_SNP;
  R["sd_ExSNP"] = SDBETA_ExSNP;
  R["cov_E_SNP"] =  COVBETA_E_SNP;
  R["cov_E_ExSNP"] = COVBETA_E_ExSNP;
  R["cov_SNP_ExSNP"] = COVBETA_SNP_ExSNP;
  return R;
}

List GxE_logit_wald_2df(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                       int beg, int end, double tol) {
  int n = Y.size();
  int r = X.ncol();

  // recopiage des matrices... en float
  /*
  MatrixXf y(n,1);
  MatrixXf x(n,r);
  for(int i = 0; i < n; i++) y(i,0) = (float) Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = (float) X(i,j);
  */
  
  // recopiage des matrices
  MatrixXd y(n,1);
  MatrixXd x(n,r);
  for(int i = 0; i < n; i++) y(i,0) = Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = X(i,j);

  // declare vectors containing result
  NumericVector BETA_E(end-beg+1);
  NumericVector BETA_SNP(end-beg+1);
  NumericVector BETA_ExSNP(end-beg+1);
  NumericVector SDBETA_E(end-beg+1);
  NumericVector SDBETA_SNP(end-beg+1);
  NumericVector SDBETA_ExSNP(end-beg+1);
  NumericVector COVBETA_E_SNP(end-beg+1);
  NumericVector COVBETA_E_ExSNP(end-beg+1);
  NumericVector COVBETA_SNP_ExSNP(end-beg+1);
  NumericVector W(end-beg+1);

  //VectorXf Beta(r);
  //Beta.setZero();
  VectorXd beta(r);
  beta.setZero();
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
	  W(i-beg) = NAN;
      continue;
    }
    // remplir les dernières colonnes de x par génotype au SNP et l'interaction(manquant -> mu)
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        x(4*ii+ss, r-2) = ((xx&3) != 3)?(xx&3):mu(i);
		x(4*ii+ss, r-1) = x(4*ii+ss, r-3)*x(4*ii+ss, r-2);
        xx >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        x(4*ii+ss, r-2) = ((xx&3) != 3)?(xx&3):mu(i);
		x(4*ii+ss, r-1) = x(4*ii+ss, r-3)*x(4*ii+ss, r-2);
        xx >>= 2;
      }
    }

	MatrixXd varbeta(r,r);
	logistic_model2<double>(y, x, beta, varbeta, tol);
	
    BETA_E(i-beg) = beta(r-3);
    BETA_SNP(i-beg) = beta(r-2);
    BETA_ExSNP(i-beg) = beta(r-1);
    SDBETA_E(i-beg) = sqrt(varbeta(r-3,r-3));
    SDBETA_SNP(i-beg) = sqrt(varbeta(r-2,r-2));
    SDBETA_ExSNP(i-beg) = sqrt(varbeta(r-1,r-1));
    COVBETA_E_SNP(i-beg) = varbeta(r-3,r-2);
    COVBETA_E_ExSNP(i-beg) = varbeta(r-3,r-1);
    COVBETA_SNP_ExSNP(i-beg) = varbeta(r-2,r-1);
	  
	VectorXd beta_int = beta.tail(2);
	MatrixXd VAR = varbeta.bottomRightCorner(2,2);
	double wald;
	wald_compute(beta_int, VAR, wald);
	W(i-beg) = wald;
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
  R["Wald"] = W;
  return R;
}


List GxE_logit_wald_3df(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                       int beg, int end, double tol) {
  int n = Y.size();
  int r = X.ncol();

  // recopiage des matrices... en float
  /*
  MatrixXf y(n,1);
  MatrixXf x(n,r);
  for(int i = 0; i < n; i++) y(i,0) = (float) Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = (float) X(i,j);
  */
  
  // recopiage des matrices
  MatrixXd y(n,1);
  MatrixXd x(n,r);
  for(int i = 0; i < n; i++) y(i,0) = Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = X(i,j);

  // declare vectors containing result
  NumericVector BETA_E(end-beg+1);
  NumericVector BETA_SNP(end-beg+1);
  NumericVector BETA_ExSNP(end-beg+1);
  NumericVector SDBETA_E(end-beg+1);
  NumericVector SDBETA_SNP(end-beg+1);
  NumericVector SDBETA_ExSNP(end-beg+1);
  NumericVector COVBETA_E_SNP(end-beg+1);
  NumericVector COVBETA_E_ExSNP(end-beg+1);
  NumericVector COVBETA_SNP_ExSNP(end-beg+1);
  NumericVector W(end-beg+1);

  //VectorXf Beta(r);
  //Beta.setZero();
  VectorXd beta(r);
  beta.setZero();
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
	  W(i-beg) = NAN;
      continue;
    }
    // remplir les dernières colonnes de x par génotype au SNP et l'interaction(manquant -> mu)
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        x(4*ii+ss, r-2) = ((xx&3) != 3)?(xx&3):mu(i);
		x(4*ii+ss, r-1) = x(4*ii+ss, r-3)*x(4*ii+ss, r-2);
        xx >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        x(4*ii+ss, r-2) = ((xx&3) != 3)?(xx&3):mu(i);
		x(4*ii+ss, r-1) = x(4*ii+ss, r-3)*x(4*ii+ss, r-2);
        xx >>= 2;
      }
    }

	MatrixXd varbeta(r,r);
	logistic_model2<double>(y, x, beta, varbeta, tol);
	
    BETA_E(i-beg) = beta(r-3);
    BETA_SNP(i-beg) = beta(r-2);
    BETA_ExSNP(i-beg) = beta(r-1);
    SDBETA_E(i-beg) = sqrt(varbeta(r-3,r-3));
    SDBETA_SNP(i-beg) = sqrt(varbeta(r-2,r-2));
    SDBETA_ExSNP(i-beg) = sqrt(varbeta(r-1,r-1));
    COVBETA_E_SNP(i-beg) = varbeta(r-3,r-2);
    COVBETA_E_ExSNP(i-beg) = varbeta(r-3,r-1);
    COVBETA_SNP_ExSNP(i-beg) = varbeta(r-2,r-1);
	  
	VectorXd beta_int = beta.tail(3);
	MatrixXd VAR = varbeta.bottomRightCorner(3,3);
	double wald;
	wald_compute(beta_int, VAR, wald);
	W(i-beg) = wald;
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
  R["Wald"] = W;
  return R;
}

RcppExport SEXP gg_GxE_logit_wald_1df(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_logit_wald_1df(pA, mu, Y, X, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_GxE_logit_wald_2df(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_logit_wald_2df(pA, mu, Y, X, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_GxE_logit_wald_3df(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_logit_wald_3df(pA, mu, Y, X, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}
