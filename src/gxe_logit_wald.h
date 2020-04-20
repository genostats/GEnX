#include <Rcpp.h>
#include "gaston/logit_model.h"
#include "gaston/matrix4.h"
#include "gaston.utils/snp_filler.h"
#include "wald.h"
#include <ctime>
#include <cmath>
#include <iostream>
#define BLOCK 20 

#define MATRIX MatrixXd
#define VECTOR VectorXd

using namespace Rcpp;
using namespace Eigen;


class GxE_logit_wald {
public:
  int n, r;
  MATRIX x, y;
  VECTOR SNP;
  double tol;
  snp_filler<double> &S;
  std::vector<double> BETA_E, BETA_SNP, BETA_ExSNP;
  std::vector<double> SDBETA_E, SDBETA_SNP, SDBETA_ExSNP;
  std::vector<double> COVBETA_E_SNP, COVBETA_E_ExSNP, COVBETA_SNP_ExSNP;
  std::vector<double> W;
  
  GxE_logit_wald(NumericVector Y, NumericMatrix X, double tol_, snp_filler<double> &S_ ) : 
                  n(Y.size()), r(X.ncol()), x(n, r), y(n,1), SNP(n), tol(tol_), S(S_) {
  
  for(int i = 0; i < n; i++) y(i,0) = Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = X(i,j);
  };
    
  void run_tests_1df() {
    VECTOR beta(r);
    MATRIX varbeta(r,r);
    beta.setZero();
      
    while( S.snp_fill(&SNP[0]) ) {
      x.col(r-2) = SNP;
      x.col(r-1) = SNP.cwiseProduct(x.col(r-3));
      
      // Model
      logistic_model2<double>(y, x, beta, varbeta, tol);
      
      BETA_E.push_back( beta(r-3) );
      BETA_SNP.push_back( beta(r-2) );
      BETA_ExSNP.push_back( beta(r-1) );
      SDBETA_E.push_back( sqrt(varbeta(r-3,r-3)) );
      SDBETA_SNP.push_back( sqrt(varbeta(r-2,r-2)) );
      SDBETA_ExSNP.push_back( sqrt(varbeta(r-1,r-1)) );
      COVBETA_E_SNP.push_back( varbeta(r-3,r-2) );
      COVBETA_E_ExSNP.push_back( varbeta(r-3,r-1) );
      COVBETA_SNP_ExSNP.push_back( varbeta(r-2,r-1) );
      W.push_back( beta(r-1)*beta(r-1)/varbeta(r-1,r-1) );
    }
    
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
    VECTOR beta(r), beta_int(2);
    MATRIX varbeta(r,r), VAR(2,2);
    double wald;
    beta.setZero();
      
    while( S.snp_fill(&SNP[0]) ) {
      x.col(r-2) = SNP;
      x.col(r-1) = SNP.cwiseProduct(x.col(r-3));
      
      // Model
      logistic_model2<double>(y, x, beta, varbeta, tol);
      
      BETA_E.push_back( beta(r-3) );
      BETA_SNP.push_back( beta(r-2) );
      BETA_ExSNP.push_back( beta(r-1) );
      SDBETA_E.push_back( sqrt(varbeta(r-3,r-3)) );
      SDBETA_SNP.push_back( sqrt(varbeta(r-2,r-2)) );
      SDBETA_ExSNP.push_back( sqrt(varbeta(r-1,r-1)) );
      COVBETA_E_SNP.push_back( varbeta(r-3,r-2) );
      COVBETA_E_ExSNP.push_back( varbeta(r-3,r-1) );
      COVBETA_SNP_ExSNP.push_back( varbeta(r-2,r-1) );
      beta_int = beta.tail(2);
      VAR = varbeta.bottomRightCorner(2,2);
      wald_compute(beta_int, VAR, wald);
      W.push_back( wald );
    }
    
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
    VECTOR beta(r), beta_int(3);
    MATRIX varbeta(r,r), VAR(3,3);
    double wald;
    beta.setZero();
      
    while( S.snp_fill(&SNP[0]) ) {
      x.col(r-2) = SNP;
      x.col(r-1) = SNP.cwiseProduct(x.col(r-3));
      
      // Model
      logistic_model2<double>(y, x, beta, varbeta, tol);
      
      BETA_E.push_back( beta(r-3) );
      BETA_SNP.push_back( beta(r-2) );
      BETA_ExSNP.push_back( beta(r-1) );
      SDBETA_E.push_back( sqrt(varbeta(r-3,r-3)) );
      SDBETA_SNP.push_back( sqrt(varbeta(r-2,r-2)) );
      SDBETA_ExSNP.push_back( sqrt(varbeta(r-1,r-1)) );
      COVBETA_E_SNP.push_back( varbeta(r-3,r-2) );
      COVBETA_E_ExSNP.push_back( varbeta(r-3,r-1) );
      COVBETA_SNP_ExSNP.push_back( varbeta(r-2,r-1) );
      beta_int = beta.tail(3);
      VAR = varbeta.bottomRightCorner(3,3);
      wald_compute(beta_int, VAR, wald);
      W.push_back( wald );
    }
    
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
