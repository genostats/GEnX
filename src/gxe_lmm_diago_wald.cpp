#include <Rcpp.h>
#include "diago2.h"
#include "matrix4.h"
#include "wald.h"
#include <cmath>

// laisser en double ça va aussi vite (plus vite ?) et ça fait vraiment
// une différence si il y a des covariables
#define scalar double
#define MATRIX MatrixXd
#define VECTOR VectorXd

//[[Rcpp::export]]
List GxE_lmm_wald(XPtr<matrix4> pA, NumericVector mu, NumericVector Y, NumericMatrix X, 
                  int p, NumericVector Sigma, NumericMatrix U, int beg, int end, double tol) {

  int n = Sigma.size();
  int r = X.ncol();
  int max_iter = 10;

  if(Y.size() != n || X.nrow() != n || U.nrow() != n || U.ncol() != n) 
    stop("Dimensions mismatch");
    
  // conversion en float si nécessaire... sinon copie
  MATRIX y0(n, 1);
  for(int i = 0; i < n; i++)
    y0(i,0) = Y[i];

  MATRIX x0(n, r);
  for(int j = 0; j < r; j++) 
    for(int i = 0; i < n; i++)
      x0(i,j) = X(i,j);

  VECTOR sigma(n);
  for(int i = 0; i < n; i++) 
    sigma[i] = Sigma[i];

  MATRIX u(n,n);
  for(int j = 0; j < n; j++) 
    for(int i = 0; i < n; i++)
      u(i,j) = U(i,j);

  MATRIX x = u.transpose() * x0;
  MATRIX y = u.transpose() * y0;

  // Vecteur SNPs and interactions
  VECTOR SNP(n);
  VECTOR INT(n);

  // declare vectors containing result
  NumericVector H2(end-beg+1);
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
  
  // variance matrix for Wald test
  VECTOR beta_int;
  MATRIX VAR;
  scalar wald;
  
  // object for likelihood maximization
  diag_likelihood<MATRIX, VECTOR, scalar> A(p, y, x, sigma);

  scalar h2 = 0;

  for(int i = beg; i <= end; i++) {
    // if(!(i%65536)) Rcout << "i = " << i << "\n";
    if( std::isnan(mu(i)) || mu(i) == 0 || mu(i) == 2 ) {
      H2(i-beg) = NAN;
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

    // remplir dernière colonne de x : récupérer SNP, multiplier par u'...
    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
		INT(4*ii+ss) = A.X(4*ii+ss,r-3) * SNP(4*ii+ss);
        x >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t x = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu(i);
		INT(4*ii+ss) = A.X(4*ii+ss,r-3) * SNP(4*ii+ss);
        x >>= 2;
      }
    }
    A.X.col(r-2) = u.transpose() * SNP;
	A.X.col(r-1) = u.transpose() * INT;

    // likelihood maximization
    h2 = (h2 > 0.9)?0.9:h2;
    A.newton_max( h2, 0, 0.99, tol, max_iter, false);
    
    // CALCUL DES BLUPS 
    VECTOR beta, omega;
    A.blup(h2, beta, omega, false, true);

// Rcout << "beta = " << beta.transpose() << "\n";
// Rcout << "v = " << A.v << "\n";
// Rcout << "XViX = " << A.XViX << "\n";
// Rcout << "XViX_i = " << A.XViX_i << "\n";
    
    if(A.d != 0) {
      H2(i-beg) = h2;
	  VAR = A.v * A.XViX_i.block(r-3,r-3,3,3);
      BETA_E(i-beg) = beta(r-3);
      BETA_SNP(i-beg) = beta(r-2);
      BETA_ExSNP(i-beg) = beta(r-1);
      SDBETA_E(i-beg) = sqrt(VAR(0,0));
      SDBETA_SNP(i-beg) = sqrt(VAR(1,1));
      SDBETA_ExSNP(i-beg) = sqrt(VAR(2,2));
      COVBETA_E_SNP(i-beg) = VAR(0,1);
      COVBETA_E_ExSNP(i-beg) = VAR(0,2);
      COVBETA_SNP_ExSNP(i-beg) = VAR(1,2);
	  
	  beta_int = beta.tail(3);
	  wald_compute(beta_int, VAR, wald);
	  WALD_3df(i-beg) = wald;
	  
	  beta_int = beta.tail(2);
	  VAR = VAR.block(1,1,2,2);
	  wald_compute(beta_int, VAR, wald);
	  WALD_2df(i-beg) = wald;
    } else {
      H2(i-beg) = NAN;
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
    }
  }

  List R;
  R["h2"] = H2;
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


RcppExport SEXP gg_GxE_lmm_wald(SEXP pASEXP, SEXP muSEXP, SEXP YSEXP, SEXP XSEXP, SEXP pSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP );
        Rcpp::traits::input_parameter< int >::type beg(begSEXP );
        Rcpp::traits::input_parameter< int >::type end(endSEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        List __result = GxE_lmm_wald(pA, mu, Y, X, p, Sigma, U, beg, end, tol);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}



