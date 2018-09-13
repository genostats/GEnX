#include <Rcpp.h>
#include "gaston/ai-reml-logit-1k-covar.h"
#include "gaston/matrix4.h"
#include "wald.h"
#include <ctime>
#include <cmath>
#define BLOCK 20 

//[[Rcpp::export]]
List GxE_logitmm_wald_1df(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                       NumericMatrix K, int beg, int end, double tol) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  Map_MatrixXd kk(as<Map<MatrixXd> >(K));

  int n = y.rows();
  int r = x.cols();
  
  /*
  // recopiage des matrices... en float
  MatrixXf y(n,1);
  MatrixXf x(n,r);
  MatrixXf kk(n,n);
  for(int i = 0; i < n; i++) y(i,0) = (float) Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = (float) X(i,j);

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++)
      kk(i,j) = (float) K(i,j);
  */

  // declare vectors containing result
  VectorXd TAU(end-beg+1);
  VectorXd BETA_E(end-beg+1);
  VectorXd BETA_SNP(end-beg+1);
  VectorXd BETA_ExSNP(end-beg+1);
  VectorXd SDBETA_E(end-beg+1);
  VectorXd SDBETA_SNP(end-beg+1);
  VectorXd SDBETA_ExSNP(end-beg+1);
  VectorXd COVBETA_E_SNP(end-beg+1);
  VectorXd COVBETA_E_ExSNP(end-beg+1);
  VectorXd COVBETA_SNP_ExSNP(end-beg+1);

  // initial values for beta, tau
  double tau = 0; 
  int niter;
  MatrixXd P(n,n);
  VectorXd omega(n);
  VectorXd beta(r);
  MatrixXd varbeta(r,r);
  
  /*
  VectorXf TAU(end-beg+1);
  VectorXf BETA_E(end-beg+1);
  VectorXf BETA_SNP(end-beg+1);
  VectorXf BETA_ExSNP(end-beg+1);
  VectorXf SDBETA_E(end-beg+1);
  VectorXf SDBETA_SNP(end-beg+1);
  VectorXf SDBETA_ExSNP(end-beg+1);
  VectorXf COVBETA_E_SNP(end-beg+1);
  VectorXf COVBETA_E_ExSNP(end-beg+1);
  VectorXf COVBETA_SNP_ExSNP(end-beg+1);

  // initial values for beta, tau
  float tau; 
  int niter;
  MatrixXf P(n,n);
  VectorXf omega(n);
  VectorXf beta(r);
  MatrixXf varbeta(r,r);
  
  for(int i = 0; i < n; i++) {
    x(i,r-2) = 0;
    x(i,r-1) = 0; }
  AIREML1_logit_f(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, false, false);
  */

  for(int i = 0; i < n; i++) {
	  x(i,r-2) = 0;
      x(i,r-1) = 0; }
  AIREML1_logit(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, false, false);

  // Rcout << min_h2 << " < h2 < " << max_h2 << "\n";
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      TAU(i-beg) = NAN;
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
    // remplir dernière colonne de x par génotype au SNP (manquant -> mu)
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

    // use last computed tau as starting point...
    if( std::isnan(tau) ) tau = 0;
    AIREML1_logit(y, x, kk, true, 0, 50, tol, false, tau, niter, P, omega, beta, varbeta, true, false);
	//AIREML1_logit_f(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, true, true);

    TAU(i-beg) = tau;
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
  R["tau"] = TAU;
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

//[[Rcpp::export]]
List GxE_logitmm_wald_2df(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                       NumericMatrix K, int beg, int end, double tol) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  Map_MatrixXd kk(as<Map<MatrixXd> >(K));

  int n = y.rows();
  int r = x.cols();
  
  /*
  // recopiage des matrices... en float
  MatrixXf y(n,1);
  MatrixXf x(n,r);
  MatrixXf kk(n,n);
  for(int i = 0; i < n; i++) y(i,0) = (float) Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = (float) X(i,j);

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++)
      kk(i,j) = (float) K(i,j);
  */

  // declare vectors containing result
  VectorXd TAU(end-beg+1);
  VectorXd BETA_E(end-beg+1);
  VectorXd BETA_SNP(end-beg+1);
  VectorXd BETA_ExSNP(end-beg+1);
  VectorXd SDBETA_E(end-beg+1);
  VectorXd SDBETA_SNP(end-beg+1);
  VectorXd SDBETA_ExSNP(end-beg+1);
  VectorXd COVBETA_E_SNP(end-beg+1);
  VectorXd COVBETA_E_ExSNP(end-beg+1);
  VectorXd COVBETA_SNP_ExSNP(end-beg+1);
  VectorXd W(end-beg+1);

  // initial values for beta, tau
  double tau = 0; 
  int niter;
  MatrixXd P(n,n);
  VectorXd omega(n);
  VectorXd beta(r);
  MatrixXd varbeta(r,r);
  
  /*
  VectorXf TAU(end-beg+1);
  VectorXf BETA_E(end-beg+1);
  VectorXf BETA_SNP(end-beg+1);
  VectorXf BETA_ExSNP(end-beg+1);
  VectorXf SDBETA_E(end-beg+1);
  VectorXf SDBETA_SNP(end-beg+1);
  VectorXf SDBETA_ExSNP(end-beg+1);
  VectorXf COVBETA_E_SNP(end-beg+1);
  VectorXf COVBETA_E_ExSNP(end-beg+1);
  VectorXf COVBETA_SNP_ExSNP(end-beg+1);
  VectorXf W(end-beg+1);

  // initial values for beta, tau
  float tau; 
  int niter;
  MatrixXf P(n,n);
  VectorXf omega(n);
  VectorXf beta(r);
  MatrixXf varbeta(r,r);
  
  for(int i = 0; i < n; i++) {
    x(i,r-2) = 0;
    x(i,r-1) = 0; }
  AIREML1_logit_f(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, false, false);
  */

  for(int i = 0; i < n; i++) {
	  x(i,r-2) = 0;
      x(i,r-1) = 0; }
  AIREML1_logit(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, false, false);

  // Rcout << min_h2 << " < h2 < " << max_h2 << "\n";
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      TAU(i-beg) = NAN;
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
    // remplir dernière colonne de x par génotype au SNP (manquant -> mu)
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

    // use last computed tau as starting point...
    if( std::isnan(tau) ) tau = 0;
    AIREML1_logit(y, x, kk, true, 0, 50, tol, false, tau, niter, P, omega, beta, varbeta, true, false);
	//AIREML1_logit_f(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, true, true);

    TAU(i-beg) = tau;
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
  R["tau"] = TAU;
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


//[[Rcpp::export]]
List GxE_logitmm_wald_3df(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                       NumericMatrix K, int beg, int end, double tol) {
  Map_MatrixXd y(as<Map<MatrixXd> >(Y));
  Map_MatrixXd x(as<Map<MatrixXd> >(X));
  Map_MatrixXd kk(as<Map<MatrixXd> >(K));

  int n = y.rows();
  int r = x.cols();
  
  /*
  // recopiage des matrices... en float
  MatrixXf y(n,1);
  MatrixXf x(n,r);
  MatrixXf kk(n,n);
  for(int i = 0; i < n; i++) y(i,0) = (float) Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = (float) X(i,j);

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++)
      kk(i,j) = (float) K(i,j);
  */

  // declare vectors containing result
  VectorXd TAU(end-beg+1);
  VectorXd BETA_E(end-beg+1);
  VectorXd BETA_SNP(end-beg+1);
  VectorXd BETA_ExSNP(end-beg+1);
  VectorXd SDBETA_E(end-beg+1);
  VectorXd SDBETA_SNP(end-beg+1);
  VectorXd SDBETA_ExSNP(end-beg+1);
  VectorXd COVBETA_E_SNP(end-beg+1);
  VectorXd COVBETA_E_ExSNP(end-beg+1);
  VectorXd COVBETA_SNP_ExSNP(end-beg+1);
  VectorXd W(end-beg+1);

  // initial values for beta, tau
  double tau = 0; 
  int niter;
  MatrixXd P(n,n);
  VectorXd omega(n);
  VectorXd beta(r);
  MatrixXd varbeta(r,r);
  
  /*
  VectorXf TAU(end-beg+1);
  VectorXf BETA_E(end-beg+1);
  VectorXf BETA_SNP(end-beg+1);
  VectorXf BETA_ExSNP(end-beg+1);
  VectorXf SDBETA_E(end-beg+1);
  VectorXf SDBETA_SNP(end-beg+1);
  VectorXf SDBETA_ExSNP(end-beg+1);
  VectorXf COVBETA_E_SNP(end-beg+1);
  VectorXf COVBETA_E_ExSNP(end-beg+1);
  VectorXf COVBETA_SNP_ExSNP(end-beg+1);
  VectorXf W(end-beg+1);

  // initial values for beta, tau
  float tau; 
  int niter;
  MatrixXf P(n,n);
  VectorXf omega(n);
  VectorXf beta(r);
  MatrixXf varbeta(r,r);
  
  for(int i = 0; i < n; i++) {
    x(i,r-2) = 0;
    x(i,r-1) = 0; }
  AIREML1_logit_f(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, false, false);
  */

  for(int i = 0; i < n; i++) {
	  x(i,r-2) = 0;
      x(i,r-1) = 0; }
  AIREML1_logit(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, false, false);

  // Rcout << min_h2 << " < h2 < " << max_h2 << "\n";
  for(int i = beg; i <= end; i++) {
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      TAU(i-beg) = NAN;
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
    // remplir dernière colonne de x par génotype au SNP (manquant -> mu)
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

    // use last computed tau as starting point...
    if( std::isnan(tau) ) tau = 0;
    AIREML1_logit(y, x, kk, true, 0, 50, tol, false, tau, niter, P, omega, beta, varbeta, true, false);
	//AIREML1_logit_f(y, x, kk, true, 1e-6, 100, tol, false, tau, niter, P, omega, beta, varbeta, true, true);

    TAU(i-beg) = tau;
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
  R["tau"] = TAU;
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



RcppExport SEXP gg_GxE_logitmm_wald_1df(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP KSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_logitmm_wald_1df(pA, mu, Y, X, K, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_GxE_logitmm_wald_2df(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP KSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_logitmm_wald_2df(pA, mu, Y, X, K, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_GxE_logitmm_wald_3df(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP KSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GxE_logitmm_wald_3df(pA, mu, Y, X, K, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

