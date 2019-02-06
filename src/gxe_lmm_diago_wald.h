#include <Rcpp.h>
#include "gaston/diago2.h"
#include "gaston.utils/snp_filler.h"
#include "wald.h"
#include <cmath>

// laisser en double ça va aussi vite (plus vite ?) et ça fait vraiment
// une différence si il y a des covariables
#define MATRIX MatrixXd
#define VECTOR VectorXd

class GxE_lmm_wald { 
public:
  int n, r, p;
  MATRIX x0, u, y0;
  VECTOR sigma, SNP;
  double tol;
  snp_filler<double> &S;
  std::vector<double> H2;
  std::vector<double> BETA_E, BETA_SNP, BETA_ExSNP;
  std::vector<double> SDBETA_E, SDBETA_SNP, SDBETA_ExSNP;
  std::vector<double> COVBETA_E_SNP, COVBETA_E_ExSNP, COVBETA_SNP_ExSNP;
  std::vector<double> W;
  MATRIX x, y;
  
  GxE_lmm_wald(NumericVector Y, NumericMatrix X, 
                  int p_, NumericVector Sigma, NumericMatrix U, double tol_,
                  snp_filler<double> &S_ ) : 
                  n(Sigma.size()), r(X.ncol()), p(p_), x0(n, r), u(n,n),
                  y0(n,1), sigma(n), SNP(n), tol(tol_), S(S_) {

    if(Y.size() != n || X.nrow() != n || U.nrow() != n || U.ncol() != n) 
      stop("Dimensions mismatch");
    
    for(int i = 0; i < n; i++)
      y0(i,0) = Y[i];

    for(int j = 0; j < r; j++) 
      for(int i = 0; i < n; i++)
        x0(i,j) = X(i,j);

    for(int i = 0; i < n; i++) 
      sigma[i] = Sigma[i];

    for(int j = 0; j < n; j++) 
      for(int i = 0; i < n; i++)
        u(i,j) = U(i,j);

    x = u.transpose() * x0;
    y = u.transpose() * y0;
  }
  
  void run_tests_1df() {
    int max_iter=10;
    double h2 = 0;
    VECTOR beta, omega;
    MATRIX VAR;
    
    // object for likelihood maximization
    diag_likelihood<MATRIX, VECTOR, double> A(p, y, x, sigma);
    
    while( S.snp_fill(&SNP[0]) ) {
      A.X.col(r-2) = u.transpose() * SNP;
	  A.X.col(r-1) = u.transpose() * SNP.cwiseProduct(x0.col(r-3));

      // likelihood maximization
      h2 = (h2 > 0.9)?0.9:h2;
      A.newton_max( h2, 0, 0.99, tol, max_iter, false);

      // CALCUL DES BLUPS 
      A.blup(h2, beta, omega, false, true);

      if(A.d != 0) {
        H2.push_back(h2);
        VAR = A.v * A.XViX_i.bottomRightCorner(3,3);
        BETA_E.push_back( beta(r-3) );
        BETA_SNP.push_back( beta(r-2) );
        BETA_ExSNP.push_back( beta(r-1) );
        SDBETA_E.push_back( sqrt(VAR(0,0)) );
        SDBETA_SNP.push_back( sqrt(VAR(1,1)) );
        SDBETA_ExSNP.push_back( sqrt(VAR(2,2)) );
        COVBETA_E_SNP.push_back( VAR(0,1) );
        COVBETA_E_ExSNP.push_back( VAR(0,2) );
        COVBETA_SNP_ExSNP.push_back( VAR(1,2) );
        W.push_back( beta(r-1)*beta(r-1)/VAR(2,2) );
      } else {
        H2.push_back( NAN );
        BETA_E.push_back( NAN );
        BETA_SNP.push_back( NAN );
        BETA_ExSNP.push_back( NAN );
        SDBETA_E.push_back( NAN );
        SDBETA_SNP.push_back( NAN );
        SDBETA_ExSNP.push_back( NAN );
        COVBETA_E_SNP.push_back( NAN );
        COVBETA_E_ExSNP.push_back( NAN );
        COVBETA_SNP_ExSNP.push_back( NAN );
        W.push_back( NAN );
      }
    }
    
    S.L["h2"] = wrap( H2 );
    S.L["beta_E"] = wrap( BETA_E );
    S.L["beta_SNP"] = wrap( BETA_SNP );
    S.L["beta_ExSNP"] = wrap( BETA_ExSNP );
    S.L["sd_E"] = wrap( SDBETA_E );
    S.L["sd_SNP"] = wrap( SDBETA_SNP );
    S.L["sd_ExSNP"] = wrap( SDBETA_ExSNP );
    S.L["cov_E_SNP"] =  wrap( COVBETA_E_SNP );
    S.L["cov_E_ExSNP"] = wrap( COVBETA_E_ExSNP );
    S.L["cov_SNP_ExSNP"] = wrap( COVBETA_SNP_ExSNP );
    S.L["Wald"] = wrap( W );
  };
  
  void run_tests_2df() {
    int max_iter=10;
    double h2 = 0;
    VECTOR beta, omega;
    MATRIX VAR;
    double wald;
    
    // object for likelihood maximization
    diag_likelihood<MATRIX, VECTOR, double> A(p, y, x, sigma);
    
    while( S.snp_fill(&SNP[0]) ) {
      A.X.col(r-2) = u.transpose() * SNP;
	  A.X.col(r-1) = u.transpose() * SNP.cwiseProduct(x0.col(r-3));

      // likelihood maximization
      h2 = (h2 > 0.9)?0.9:h2;
      A.newton_max( h2, 0, 0.99, tol, max_iter, false);

      // CALCUL DES BLUPS 
      A.blup(h2, beta, omega, false, true);

      if(A.d != 0) {
        H2.push_back(h2);
        VAR = A.v * A.XViX_i.bottomRightCorner(3,3);
        BETA_E.push_back( beta(r-3) );
        BETA_SNP.push_back( beta(r-2) );
        BETA_ExSNP.push_back( beta(r-1) );
        SDBETA_E.push_back( sqrt(VAR(0,0)) );
        SDBETA_SNP.push_back( sqrt(VAR(1,1)) );
        SDBETA_ExSNP.push_back( sqrt(VAR(2,2)) );
        COVBETA_E_SNP.push_back( VAR(0,1) );
        COVBETA_E_ExSNP.push_back( VAR(0,2) );
        COVBETA_SNP_ExSNP.push_back( VAR(1,2) );
        
        VECTOR beta_int = beta.tail(2);
        VAR = A.v * A.XViX_i.bottomRightCorner(2,2);
        wald_compute(beta_int, VAR, wald);
	    W.push_back(wald);
      } else {
        H2.push_back( NAN );
        BETA_E.push_back( NAN );
        BETA_SNP.push_back( NAN );
        BETA_ExSNP.push_back( NAN );
        SDBETA_E.push_back( NAN );
        SDBETA_SNP.push_back( NAN );
        SDBETA_ExSNP.push_back( NAN );
        COVBETA_E_SNP.push_back( NAN );
        COVBETA_E_ExSNP.push_back( NAN );
        COVBETA_SNP_ExSNP.push_back( NAN );
        W.push_back( NAN );
      }
    }
    
    S.L["h2"] = wrap( H2 );
    S.L["beta_E"] = wrap( BETA_E );
    S.L["beta_SNP"] = wrap( BETA_SNP );
    S.L["beta_ExSNP"] = wrap( BETA_ExSNP );
    S.L["sd_E"] = wrap( SDBETA_E );
    S.L["sd_SNP"] = wrap( SDBETA_SNP );
    S.L["sd_ExSNP"] = wrap( SDBETA_ExSNP );
    S.L["cov_E_SNP"] =  wrap( COVBETA_E_SNP );
    S.L["cov_E_ExSNP"] = wrap( COVBETA_E_ExSNP );
    S.L["cov_SNP_ExSNP"] = wrap( COVBETA_SNP_ExSNP );
    S.L["Wald"] = wrap( W );
  };
  
  void run_tests_3df() {
    int max_iter=10;
    double h2 = 0;
    VECTOR beta, omega;
    MATRIX VAR;
    double wald;
    
    // object for likelihood maximization
    diag_likelihood<MATRIX, VECTOR, double> A(p, y, x, sigma);
    
    while( S.snp_fill(&SNP[0]) ) {
      A.X.col(r-2) = u.transpose() * SNP;
	  A.X.col(r-1) = u.transpose() * SNP.cwiseProduct(x0.col(r-3));

      // likelihood maximization
      h2 = (h2 > 0.9)?0.9:h2;
      A.newton_max( h2, 0, 0.99, tol, max_iter, false);

      // CALCUL DES BLUPS 
      A.blup(h2, beta, omega, false, true);

      if(A.d != 0) {
        H2.push_back(h2);
        VAR = A.v * A.XViX_i.bottomRightCorner(3,3);
        BETA_E.push_back( beta(r-3) );
        BETA_SNP.push_back( beta(r-2) );
        BETA_ExSNP.push_back( beta(r-1) );
        SDBETA_E.push_back( sqrt(VAR(0,0)) );
        SDBETA_SNP.push_back( sqrt(VAR(1,1)) );
        SDBETA_ExSNP.push_back( sqrt(VAR(2,2)) );
        COVBETA_E_SNP.push_back( VAR(0,1) );
        COVBETA_E_ExSNP.push_back( VAR(0,2) );
        COVBETA_SNP_ExSNP.push_back( VAR(1,2) );
        
        VECTOR beta_int = beta.tail(3);
        wald_compute(beta_int, VAR, wald);
        W.push_back(wald);
      } else {
        H2.push_back( NAN );
        BETA_E.push_back( NAN );
        BETA_SNP.push_back( NAN );
        BETA_ExSNP.push_back( NAN );
        SDBETA_E.push_back( NAN );
        SDBETA_SNP.push_back( NAN );
        SDBETA_ExSNP.push_back( NAN );
        COVBETA_E_SNP.push_back( NAN );
        COVBETA_E_ExSNP.push_back( NAN );
        COVBETA_SNP_ExSNP.push_back( NAN );
        W.push_back( NAN );
      }
    }
    
    S.L["h2"] = wrap( H2 );
    S.L["beta_E"] = wrap( BETA_E );
    S.L["beta_SNP"] = wrap( BETA_SNP );
    S.L["beta_ExSNP"] = wrap( BETA_ExSNP );
    S.L["sd_E"] = wrap( SDBETA_E );
    S.L["sd_SNP"] = wrap( SDBETA_SNP );
    S.L["sd_ExSNP"] = wrap( SDBETA_ExSNP );
    S.L["cov_E_SNP"] =  wrap( COVBETA_E_SNP );
    S.L["cov_E_ExSNP"] = wrap( COVBETA_E_ExSNP );
    S.L["cov_SNP_ExSNP"] = wrap( COVBETA_SNP_ExSNP );
    S.L["Wald"] = wrap( W );
  };
};


