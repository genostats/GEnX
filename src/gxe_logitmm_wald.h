#include <Rcpp.h>
#include "gaston/ai-reml-logit-1k-covar.h"
#include "gaston/matrix4.h"
#include "gaston.utils/snp_filler.h"
#include "wald.h"
#include <ctime>
#include <cmath>
#define BLOCK 20 

#define MATRIX MatrixXd
#define VECTOR VectorXd

class GxE_logitmm_wald { 
public:
  int n, r;
  MATRIX x0, y0, kk;
  VECTOR SNP;
  double tol;
  snp_filler<double> &S;
  std::vector<double> TAU;
  std::vector<double> BETA_E, BETA_SNP, BETA_ExSNP;
  std::vector<double> SDBETA_E, SDBETA_SNP, SDBETA_ExSNP;
  std::vector<double> COVBETA_E_SNP, COVBETA_E_ExSNP, COVBETA_SNP_ExSNP;
  std::vector<double> W;  

  GxE_logitmm_wald(NumericVector Y, NumericMatrix X, NumericMatrix K, double tol_,
                  snp_filler<double> &S_ ) : 
                  n(Y.size()), r(X.ncol()), x0(n, r), y0(n,1), kk(n,n), SNP(n), tol(tol_), S(S_) {

    if(Y.size() != n || X.nrow() != n || K.nrow() != n || K.ncol() != n) 
      stop("Dimensions mismatch");
    
    for(int i = 0; i < n; i++)
      y0(i,0) = Y[i];

    for(int j = 0; j < r; j++) 
      for(int i = 0; i < n; i++)
        x0(i,j) = X(i,j);

    for(int i = 0; i < n; i++) 
      for(int j = 0; j < n; j++)
        kk(i,j) = K(i,j);
  }
  
  void run_tests_1df() {
    int niter, max_iter=50;
    double tau = 0;
    VECTOR beta0(r-2), beta(r), omega(n);
    MATRIX P(n,n), xx(n,r-2), VAR0(r-2,r-2), VAR(r,r);
    
    xx=x0.leftCols(r-2);
    
    // object for likelihood maximization
    AIREML1_logit(y0, xx, kk, true, 1e-6, 10, 1e-5, false, tau, niter, P, omega, beta0, VAR0, false, false, false);
    
    while( S.snp_fill(&SNP[0]) ) {      
      x0.col(r-2) = SNP;
      x0.col(r-1) = SNP.cwiseProduct(x0.col(r-3));
    
      if( std::isnan(tau) ) tau = 0;
      AIREML1_logit(y0, x0, kk, true, 1e-6, max_iter, tol, false, tau, niter, P, omega, beta, VAR, true, false, false);

      TAU.push_back(tau);
      BETA_E.push_back( beta(r-3) );
      BETA_SNP.push_back( beta(r-2) );
      BETA_ExSNP.push_back( beta(r-1) );
      SDBETA_E.push_back( sqrt(VAR(r-3,r-3)) );
      SDBETA_SNP.push_back( sqrt(VAR(r-2,r-2)) );
      SDBETA_ExSNP.push_back( sqrt(VAR(r-1,r-1)) );
      COVBETA_E_SNP.push_back( VAR(r-3,r-2) );
      COVBETA_E_ExSNP.push_back( VAR(r-3,r-1) );
      COVBETA_SNP_ExSNP.push_back( VAR(r-2,r-1) );
      W.push_back( beta(r-1)*beta(r-1)/VAR(r-1,r-1) );
    }
    
    S.L["tau"] = wrap( TAU );
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
    int niter, max_iter=50;
    double tau = 0;
    VECTOR beta0(r-2), beta(r), omega(n);
    MATRIX P(n,n), xx(n,r-2), VAR0(r-2,r-2), VAR(r,r);
    
    xx=x0.leftCols(r-2);
    
    // object for likelihood maximization
    AIREML1_logit(y0, xx, kk, true, 1e-6, 10, 1e-5, false, tau, niter, P, omega, beta0, VAR0, false, false, false);
    
    while( S.snp_fill(&SNP[0]) ) {      
      x0.col(r-2) = SNP;
      x0.col(r-1) = SNP.cwiseProduct(x0.col(r-3));
    
      if( std::isnan(tau) ) tau = 0;
      AIREML1_logit(y0, x0, kk, true, 1e-6, max_iter, tol, false, tau, niter, P, omega, beta, VAR, true, false, false);

      TAU.push_back(tau);
      BETA_E.push_back( beta(r-3) );
      BETA_SNP.push_back( beta(r-2) );
      BETA_ExSNP.push_back( beta(r-1) );
      SDBETA_E.push_back( sqrt(VAR(r-3,r-3)) );
      SDBETA_SNP.push_back( sqrt(VAR(r-2,r-2)) );
      SDBETA_ExSNP.push_back( sqrt(VAR(r-1,r-1)) );
      COVBETA_E_SNP.push_back( VAR(r-3,r-2) );
      COVBETA_E_ExSNP.push_back( VAR(r-3,r-1) );
      COVBETA_SNP_ExSNP.push_back( VAR(r-2,r-1) );
      VectorXd beta_int = beta.tail(2);
      MatrixXd v = VAR.bottomRightCorner(2,2);
      double wald;
      wald_compute(beta_int, v, wald);
      W.push_back( wald );
    }
    
    S.L["tau"] = wrap( TAU );
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
    int niter, max_iter=50;
    double tau = 0;
    VECTOR beta0(r-2), beta(r), omega(n);
    MATRIX P(n,n), xx(n,r-2), VAR0(r-2,r-2), VAR(r,r);
    
    xx=x0.leftCols(r-2);
    
    // object for likelihood maximization
    AIREML1_logit(y0, xx, kk, true, 1e-6, 10, 1e-5, false, tau, niter, P, omega, beta0, VAR0, false, false, false);
    
    while( S.snp_fill(&SNP[0]) ) {      
      x0.col(r-2) = SNP;
      x0.col(r-1) = SNP.cwiseProduct(x0.col(r-3));
    
      if( std::isnan(tau) ) tau = 0;
      AIREML1_logit(y0, x0, kk, true, 1e-6, max_iter, tol, false, tau, niter, P, omega, beta, VAR, true, false, false);

      TAU.push_back(tau);
      BETA_E.push_back( beta(r-3) );
      BETA_SNP.push_back( beta(r-2) );
      BETA_ExSNP.push_back( beta(r-1) );
      SDBETA_E.push_back( sqrt(VAR(r-3,r-3)) );
      SDBETA_SNP.push_back( sqrt(VAR(r-2,r-2)) );
      SDBETA_ExSNP.push_back( sqrt(VAR(r-1,r-1)) );
      COVBETA_E_SNP.push_back( VAR(r-3,r-2) );
      COVBETA_E_ExSNP.push_back( VAR(r-3,r-1) );
      COVBETA_SNP_ExSNP.push_back( VAR(r-2,r-1) );
      VectorXd beta_int = beta.tail(3);
      MatrixXd v = VAR.bottomRightCorner(3,3);
      double wald;
      wald_compute(beta_int, v, wald);
      W.push_back( wald );
    }
    
    S.L["tau"] = wrap( TAU );
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
