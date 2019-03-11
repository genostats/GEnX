#include <Rcpp.h>
#include <RcppEigen.h>
#include "gaston/matrix4.h"
#include "gaston/matrix-varia.h"
#include "gaston.utils/snp_filler.h"
#include "wald.h"
#include <ctime>
#include <cmath>

#define MATRIX MatrixXd
#define VECTOR VectorXd

using namespace Rcpp;
using namespace Eigen;


class GxE_lm_quanti {
public:
  int n, r;
  MATRIX x, y;
  VECTOR SNP;
  snp_filler<double> &S;
  std::vector<double> BETA_E, BETA_SNP, BETA_ExSNP;
  std::vector<double> SDBETA_E, SDBETA_SNP, SDBETA_ExSNP;
  std::vector<double> COVBETA_E_SNP, COVBETA_E_ExSNP, COVBETA_SNP_ExSNP;
  std::vector<double> WW;
  MATRIX xtx, xtx_i, xty, xtx_ixty;
  
  GxE_lm_quanti(NumericVector Y, NumericMatrix X, snp_filler<double> &S_ ) : 
                  n(Y.size()), r(X.ncol()), x(n, r), y(n,1), SNP(n), S(S_),
                  xtx(r-3,r-3), xtx_i(r-3,r-3), xty(r-3,1), xtx_ixty(r-3,1) {

    for(int i = 0; i < n; i++) y(i,0) = Y[i];

    for(int i = 0; i < n; i++) 
      for(int j = 0; j < r; j++)
        x(i,j) = X(i,j);

    // Fixed matrix
    xtx = x.topLeftCorner(n,r-3).transpose() * x.topLeftCorner(n,r-3);
    xty = x.topLeftCorner(n,r-3).transpose() * y;
    xtx_i = xtx.llt().solve( MatrixXd::Identity(r-3,r-3) );
    xtx_ixty = xtx_i * xty;
  };
    
  void run_tests_1df() {
    // Matrix used in estimation computations
    MATRIX Ity(3,1);
    MATRIX xtI(r-3,3);
    MATRIX xtx_ixtI(r-3,3);
    MATRIX W(3,3);
    MATRIX W_i(3,3);
    MATRIX VAR(3,3);
    // declare vectors containing result
    VECTOR beta(r-3), gamma(3);
    double s2;
    
    while( S.snp_fill(&SNP[0]) ) {
      x.col(r-2) = SNP;
      x.col(r-1) = SNP.cwiseProduct(x.col(r-3));
   
      // Matrix depending on SNP
      Ity = x.bottomRightCorner(n,3).transpose() * y;
      xtI = x.topLeftCorner(n,r-3).transpose() * x.bottomRightCorner(n,3);
      xtx_ixtI = xtx_i * xtI;
      W = x.bottomRightCorner(n,3).transpose() * x.bottomRightCorner(n,3) - xtI.transpose() * ( xtx_i * xtI );
      W_i = W.llt().solve( MatrixXd::Identity(3,3) );
      beta = xtx_ixty +  xtx_ixtI * ( W_i * ( xtx_ixtI.transpose() * xty ) ) - xtx_ixtI * ( W_i * Ity );
      gamma = -W_i * ( xtx_ixtI.transpose() * xty ) + W_i * Ity;
      s2 = (y - x.topLeftCorner(n,r-3)*beta - x.bottomRightCorner(n,3)*gamma).squaredNorm()/(n-r-3);
      VAR = s2 * W_i;
      
      BETA_E.push_back( gamma(0) );
      BETA_SNP.push_back( gamma(1) );
      BETA_ExSNP.push_back( gamma(2) );
      SDBETA_E.push_back( sqrt(VAR(0,0)) );
      SDBETA_SNP.push_back( sqrt(VAR(1,1)) );
      SDBETA_ExSNP.push_back( sqrt(VAR(2,2)) );
      COVBETA_E_SNP.push_back( VAR(0,1) );
      COVBETA_E_ExSNP.push_back( VAR(0,2) );
      COVBETA_SNP_ExSNP.push_back( VAR(1,2) );
      WW.push_back( gamma(2)*gamma(2)/VAR(2,2) );
    };
    S.L["beta_E"] = wrap( BETA_E );
    S.L["beta_SNP"] = wrap( BETA_SNP );
    S.L["beta_ExSNP"] = wrap( BETA_ExSNP );
    S.L["sd_E"] = wrap( SDBETA_E );
    S.L["sd_SNP"] = wrap( SDBETA_SNP );
    S.L["sd_ExSNP"] = wrap( SDBETA_ExSNP );
    S.L["cov_E_SNP"] =  wrap( COVBETA_E_SNP );
    S.L["cov_E_ExSNP"] = wrap( COVBETA_E_ExSNP );
    S.L["cov_SNP_ExSNP"] = wrap( COVBETA_SNP_ExSNP );
    S.L["Wald"] = wrap( WW );
  };
  
  void run_tests_2df() { 
    // Matrix used in estimation computations
    MATRIX Ity(3,1);
    MATRIX xtI(r-3,3);
    MATRIX xtx_ixtI(r-3,3);
    MATRIX W(3,3);
    MATRIX W_i(3,3);
    MATRIX VAR(3,3), VAR2(2,2);
    // declare vectors containing result
    VECTOR beta(r-3), gamma(3), gamma2(2);
    double s2, wald;
    
    while( S.snp_fill(&SNP[0]) ) {
      x.col(r-2) = SNP;
      x.col(r-1) = SNP.cwiseProduct(x.col(r-3));
   
      // Matrix depending on SNP
      Ity = x.bottomRightCorner(n,3).transpose() * y;
      xtI = x.topLeftCorner(n,r-3).transpose() * x.bottomRightCorner(n,3);
      xtx_ixtI = xtx_i * xtI;
      W = x.bottomRightCorner(n,3).transpose() * x.bottomRightCorner(n,3) - xtI.transpose() * ( xtx_i * xtI );
      W_i = W.llt().solve( MatrixXd::Identity(3,3) );
      beta = xtx_ixty +  xtx_ixtI * ( W_i * ( xtx_ixtI.transpose() * xty ) ) - xtx_ixtI * ( W_i * Ity );
      gamma = -W_i * ( xtx_ixtI.transpose() * xty ) + W_i * Ity;
      s2 = (y - x.topLeftCorner(n,r-3)*beta - x.bottomRightCorner(n,3)*gamma).squaredNorm()/(n-r-3);
      VAR = s2 * W_i;
      
      BETA_E.push_back( gamma(0) );
      BETA_SNP.push_back( gamma(1) );
      BETA_ExSNP.push_back( gamma(2) );
      SDBETA_E.push_back( sqrt(VAR(0,0)) );
      SDBETA_SNP.push_back( sqrt(VAR(1,1)) );
      SDBETA_ExSNP.push_back( sqrt(VAR(2,2)) );
      COVBETA_E_SNP.push_back( VAR(0,1) );
      COVBETA_E_ExSNP.push_back( VAR(0,2) );
      COVBETA_SNP_ExSNP.push_back( VAR(1,2) );
      gamma2 = gamma.tail(2);
      VAR2 = VAR.bottomRightCorner(2,2);
      wald_compute(gamma2, VAR2, wald);
      WW.push_back( wald );
    };
    S.L["beta_E"] = wrap( BETA_E );
    S.L["beta_SNP"] = wrap( BETA_SNP );
    S.L["beta_ExSNP"] = wrap( BETA_ExSNP );
    S.L["sd_E"] = wrap( SDBETA_E );
    S.L["sd_SNP"] = wrap( SDBETA_SNP );
    S.L["sd_ExSNP"] = wrap( SDBETA_ExSNP );
    S.L["cov_E_SNP"] =  wrap( COVBETA_E_SNP );
    S.L["cov_E_ExSNP"] = wrap( COVBETA_E_ExSNP );
    S.L["cov_SNP_ExSNP"] = wrap( COVBETA_SNP_ExSNP );
    S.L["Wald"] = wrap( WW );
  };
  
  void run_tests_3df() {
    // Matrix used in estimation computations
    MATRIX Ity(3,1);
    MATRIX xtI(r-3,3);
    MATRIX xtx_ixtI(r-3,3);
    MATRIX W(3,3);
    MATRIX W_i(3,3);
    MATRIX VAR(3,3);
    // declare vectors containing result
    VECTOR beta(r-3), gamma(3);
    double s2, wald;
    
    while( S.snp_fill(&SNP[0]) ) {
      x.col(r-2) = SNP;
      x.col(r-1) = SNP.cwiseProduct(x.col(r-3));
   
      // Matrix depending on SNP
      Ity = x.bottomRightCorner(n,3).transpose() * y;
      xtI = x.topLeftCorner(n,r-3).transpose() * x.bottomRightCorner(n,3);
      xtx_ixtI = xtx_i * xtI;
      W = x.bottomRightCorner(n,3).transpose() * x.bottomRightCorner(n,3) - xtI.transpose() * ( xtx_i * xtI );
      W_i = W.llt().solve( MatrixXd::Identity(3,3) );
      beta = xtx_ixty +  xtx_ixtI * ( W_i * ( xtx_ixtI.transpose() * xty ) ) - xtx_ixtI * ( W_i * Ity );
      gamma = -W_i * ( xtx_ixtI.transpose() * xty ) + W_i * Ity;
      s2 = (y - x.topLeftCorner(n,r-3)*beta - x.bottomRightCorner(n,3)*gamma).squaredNorm()/(n-r-3);
      VAR = s2 * W_i;
      
      BETA_E.push_back( gamma(0) );
      BETA_SNP.push_back( gamma(1) );
      BETA_ExSNP.push_back( gamma(2) );
      SDBETA_E.push_back( sqrt(VAR(0,0)) );
      SDBETA_SNP.push_back( sqrt(VAR(1,1)) );
      SDBETA_ExSNP.push_back( sqrt(VAR(2,2)) );
      COVBETA_E_SNP.push_back( VAR(0,1) );
      COVBETA_E_ExSNP.push_back( VAR(0,2) );
      COVBETA_SNP_ExSNP.push_back( VAR(1,2) );
      wald_compute(gamma, VAR, wald);
      WW.push_back( wald );
    };
    S.L["beta_E"] = wrap( BETA_E );
    S.L["beta_SNP"] = wrap( BETA_SNP );
    S.L["beta_ExSNP"] = wrap( BETA_ExSNP );
    S.L["sd_E"] = wrap( SDBETA_E );
    S.L["sd_SNP"] = wrap( SDBETA_SNP );
    S.L["sd_ExSNP"] = wrap( SDBETA_ExSNP );
    S.L["cov_E_SNP"] =  wrap( COVBETA_E_SNP );
    S.L["cov_E_ExSNP"] = wrap( COVBETA_E_ExSNP );
    S.L["cov_SNP_ExSNP"] = wrap( COVBETA_SNP_ExSNP );
    S.L["Wald"] = wrap( WW );
  };
};


