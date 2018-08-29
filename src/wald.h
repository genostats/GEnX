#ifndef WALD
#define WALD
#include <RcppEigen.h>
#include <cmath>
#include <iostream>
#include "gaston/matrix-varia.h"

using namespace Rcpp;
using namespace Eigen;

inline void wald_compute(Eigen::VectorXd & beta, Eigen::MatrixXd & VAR, double & wald) {
  int r(VAR.rows());
  MatrixXd VAR_i(r,r);
  double d, log_d;

  sym_inverse(VAR, VAR_i, log_d, d, 1e-5);
  wald = beta.transpose() * ( VAR_i * beta );
}

#endif
