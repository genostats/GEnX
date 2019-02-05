#include <Rcpp.h>
#include "gaston/diago2_full.h"
#include "gaston/diago2_full_nocovar.h"
#include "gaston.utils/snp_filler.h"

#define MATRIX MatrixXd
#define VECTOR VectorXd

class GxE_lmm_lrt { 
public:
  int n, r, p;
  MATRIX x0, u, y0;
  VECTOR sigma, SNP;
  double tol;
  snp_filler<double> &S;
  std::vector<double> H2, LRT;
  MATRIX x, y;

    
  GxE_lmm_lrt(NumericVector Y, NumericMatrix X, 
                  int p_, NumericVector Sigma, NumericMatrix U, double tol_,
                  snp_filler<double> &S_ ) : 
                  n(Sigma.size()), r(X.ncol()), p(p_), x0(n, r), u(n,n),
                  y0(n,1), sigma(n), SNP(n), tol(tol_), S(S_) {
                      
    if(Y.size() != n || X.nrow() != n || U.nrow() != n || U.ncol() != n)
      stop("Dimensions mismatch");

    // conversion en float...
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
    double h2 = 0, h20 = 0;
    
    // object for likelihood maximization, model null
    MATRIX x1 = x.leftCols(r-1);
    diag_full_likelihood<MATRIX, VECTOR, double> A0(p, y, x1, sigma);

    // object for likelihood maximization
    diag_full_likelihood<MATRIX, VECTOR, double> A(p, y, x, sigma);

  
    while( S.snp_fill(&SNP[0]) ) {
      A.X.col(r-2) = u.transpose() * SNP;
      A0.X.col(r-2) = A.X.col(r-2);
      A.X.col(r-1) = u.transpose() * SNP.cwiseProduct(x0.col(r-3));;

      // model null
      h20 = (h20 > 0.9)?0.9:h20;
      A0.newton_max(h20, 0, 0.99, tol, max_iter, false);

      // likelihood maximization
      h2 = (h2 > 0.9)?0.9:h2;
      A.newton_max( h2, 0, 0.99, tol, max_iter, false);

      H2.push_back(h2);
      LRT.push_back( 2*(A.likelihood() - A0.likelihood()) );
    }

    S.L["h2"] = wrap(H2);
    S.L["LRT"] = wrap(LRT);
  };
  

  void run_tests_2df() {
    int max_iter=10;
    double h2 = 0;

    // on commence par le modèle null
    MATRIX x1 = x.leftCols(r-2);
   
    // object for likelihood maximization
    diag_full_likelihood<MATRIX, VECTOR, double> A0(p, y, x1, sigma);
    A0.newton_max(h2, 0, 0.99, tol, max_iter, false);

    // object for likelihood maximization
    diag_full_likelihood<MATRIX, VECTOR, double> A(p, y, x, sigma);

    while( S.snp_fill(&SNP[0]) ) {
      A.X.col(r-2) = u.transpose() * SNP;
	  A.X.col(r-1) = u.transpose() * SNP.cwiseProduct(x0.col(r-3));

      // likelihood maximization
      h2 = (h2 > 0.9)?0.9:h2;
      A.newton_max( h2, 0, 0.99, tol, max_iter, false);

      H2.push_back(h2);
      LRT.push_back( 2*(A.likelihood() - A0.likelihood()) );
    }

    S.L["h2"] = wrap(H2);
    S.L["LRT"] = wrap(LRT);
  }
  
  void run_tests_3df() { 
    int max_iter=10;
    double h2 = 0;
    double likelihood0;

    if(r == 3) { 
      // object for likelihood maximization
      diag_full_likelihood_nocovar<MATRIX, VECTOR, double> A(p, y, sigma);
      A.newton_max(h2, 0, 0.99, tol, max_iter, false);
      likelihood0 = A.likelihood();
    } else { 
      MATRIX x1 = x.leftCols(r-3);
      // object for likelihood maximization
      diag_full_likelihood<MATRIX, VECTOR, double> A(p, y, x1, sigma);
      A.newton_max(h2, 0, 0.99, tol, max_iter, false);
      likelihood0 = A.likelihood();
    }
    
    // object for likelihood maximization
    diag_full_likelihood<MATRIX, VECTOR, double> A(p, y, x, sigma);

    while( S.snp_fill(&SNP[0]) ) {
      A.X.col(r-2) = u.transpose() * SNP;
	  A.X.col(r-1) = u.transpose() * SNP.cwiseProduct(x0.col(r-3));

      // likelihood maximization
      h2 = (h2 > 0.9)?0.9:h2;
      A.newton_max( h2, 0, 0.99, tol, max_iter, false);

      H2.push_back(h2);
      LRT.push_back( 2*(A.likelihood() - likelihood0) );
    }

    S.L["h2"] = wrap(H2);
    S.L["LRT"] = wrap(LRT);
  }
};

