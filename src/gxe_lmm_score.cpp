#include <Rcpp.h>
#include "gaston/matrix-varia.h"
#include "gaston/matrix4.h"
#include "score.h"

// float version
//[[Rcpp::export]]
List GxE_lmm_score_1df(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector mu, int beg, int end) {
  int n = PY.size();
  if(P.nrow() != n || P.ncol() != n) 
    stop("Dimensions mismatch\n");

  VectorXf Py(n);
  for(int i = 0; i < n; i++) Py(i) = PY[i];

  MatrixXf PP(n,n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) 
      PP(i,j) = P(i,j);
     

  int r = end-beg+1;

  VectorXf SNP(n);
  NumericVector s(r);
  double t, v;  
  
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      s(i-beg) = NAN;
      continue;
    }
    // récupérer SNP
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
        x >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
        x >>= 2;
      }
    }
    
    v = (PP.selfadjointView<Lower>()*SNP).dot(SNP);
    t = SNP.dot(Py);
    s(i-beg) = t*t/v;
  }
  
  List S;
  S["score"] = s;

  return S;
}

List GxE_lmm_score_2df(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector mu, NumericVector E, int beg, int end) {
  int n = PY.size();
  if(P.nrow() != n || P.ncol() != n) 
    stop("Dimensions mismatch\n");

  VectorXd Py(n);
  for(int i = 0; i < n; i++) Py(i) = PY[i];

  MatrixXd PP(n,n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) 
      PP(i,j) = P(i,j);
     

  int r = end-beg+1;

  MatrixXd X(n,2);
  NumericVector s(r);
  MatrixXd v(2,2);
  VectorXd t(2);  
  double score;
  
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      s(i-beg) = NAN;
      continue;
    }
    // récupérer SNP
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        X(4*ii+ss,0) = ((x&3) != 3)?(x&3):mu(i);
		X(4*ii+ss,1) = E(4*ii+ss)*X(4*ii+ss,0);
        x >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        X(4*ii+ss,0) = ((x&3) != 3)?(x&3):mu(i);
		X(4*ii+ss,1) = E(4*ii+ss)*X(4*ii+ss,0);
        x >>= 2;
      }
    }
    
    v = ( X.transpose() * PP.selfadjointView<Lower>() ) * X;
    t = X.transpose() * Py;
	score_compute(t, v, score);
    s(i-beg) = score;
  }
  
  List S;
  S["score"] = s;

  return S;
}

List GxE_lmm_score_3df(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector mu, NumericVector E, int beg, int end) {
  int n = PY.size();
  if(P.nrow() != n || P.ncol() != n) 
    stop("Dimensions mismatch\n");

  VectorXd Py(n);
  for(int i = 0; i < n; i++) Py(i) = PY[i];

  MatrixXd PP(n,n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) 
      PP(i,j) = P(i,j);
     

  int r = end-beg+1;

  MatrixXd X(n,3);
  for(int i = 0; i < n; i++) X(i,0) = E[i];
  
  NumericVector s(r);
  MatrixXd v(2,2);
  VectorXd t(2);  
  double score;
  
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      s(i-beg) = NAN;
      continue;
    }
    // récupérer SNP
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        X(4*ii+ss,1) = ((x&3) != 3)?(x&3):mu(i);
		X(4*ii+ss,2) = X(4*ii+ss,0)*X(4*ii+ss,1);
        x >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        X(4*ii+ss,1) = ((x&3) != 3)?(x&3):mu(i);
		X(4*ii+ss,2) = X(4*ii+ss,0)*X(4*ii+ss,1);
        x >>= 2;
      }
    }
    
    v = ( X.transpose() * PP.selfadjointView<Lower>() ) * X;
    t = X.transpose() * Py;
	score_compute(t, v, score);
    s(i-beg) = score;
  }
  
  List S;
  S["score"] = s;

  return S;
}


RcppExport SEXP gg_GxE_lmm_score_1df(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_lmm_score_1df(pA, PY, P, mu, beg, end));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_GxE_lmm_score_2df(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP ESEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_lmm_score_2df(pA, PY, P, mu, E, beg, end));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_GxE_lmm_score_3df(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP ESEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_lmm_score_3df(pA, PY, P, mu, E, beg, end));
    return rcpp_result_gen;
END_RCPP
}
