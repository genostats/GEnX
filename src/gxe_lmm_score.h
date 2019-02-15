#include <Rcpp.h>
#include <RcppEigen.h>
#include "gaston/matrix-varia.h"
#include "gaston.utils/snp_filler.h"
#include "score.h"


template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
class GxE_lmm_score { 
public:
  int n;
  VECTOR<scalar_t> Py;
  MATRIX<scalar_t> PP;
  VECTOR<scalar_t> e, SNP;
  snp_filler<scalar_t> & S;

  GxE_lmm_score(NumericVector PY, NumericMatrix P, NumericVector E, snp_filler<double> &S_ ) : 
                  n(PY.size()), Py(n), PP(n,n), e(n), SNP(n), S(S_) {

    if(P.nrow() != n || P.ncol() != n) 
      stop("Dimensions mismatch\n");

    for(int i = 0; i < n; i++) Py(i) = PY[i];
    
    for(int i = 0; i < n; i++) e(i) = E[i];

    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
        PP(i,j) = P(i,j);
  }

//   void run_tests_1df() {
//     std::vector<double> s;
//     double t, v;  
//   
//     while( S.snp_fill(&SNP[0]) ) {
//       if( S.current_snp_monomorphic() ) {
//         s.push_back(NAN);
//         continue;
//       }
//     
//       VectorXd INT=SNP.cwiseProduct(e);
//       
//       v = (PP.template selfadjointView<Eigen::Lower>() *INT).dot(INT);
//       t = INT.dot(Py);
//       s.push_back( t*t/v );
//     }
//     S.L["score"] = wrap(s);
//   };

  void run_tests_2df() {
    std::vector<double> s;
    MatrixXd X(n,2);
    double score;
  
    while( S.snp_fill(&SNP[0]) ) {
      if( S.current_snp_monomorphic() ) {
        s.push_back(NAN);
        continue;
      }

      X.col(0) = SNP;
      X.col(1) = SNP.cwiseProduct(e);
    
      MatrixXd v = ( X.transpose() * PP.template selfadjointView<Eigen::Lower>() ) * X;
      VectorXd t = X.transpose() * Py;
      score_compute(t, v, score);

      s.push_back( score );
    }
    S.L["score"] = wrap(s);
  };
  
  void run_tests_3df() {
    std::vector<double> s;
    MatrixXd X(n,3);
    double score;
  
    while( S.snp_fill(&SNP[0]) ) {
      if( S.current_snp_monomorphic() ) {
        s.push_back(NAN);
        continue;
      }

      X.col(0) = e;
      X.col(1) = SNP;
      X.col(2) = SNP.cwiseProduct(e);
    
      MatrixXd v = ( X.transpose() * PP.template selfadjointView<Eigen::Lower>() ) * X;
      VectorXd t = X.transpose() * Py;
      score_compute(t, v, score);

      s.push_back( score );
    }
    S.L["score"] = wrap(s);
  };
};

