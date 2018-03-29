#ifndef SCORE
#define SCORE
#include <RcppEigen.h>
#include <cmath>
#include <iostream>
#include "matrix-varia.h"

using namespace Rcpp;
using namespace Eigen;

inline void score_compute(Eigen::VectorXd & beta, Eigen::MatrixXd & VAR, double & score) {
  int r(VAR.rows());
  MatrixXd VAR_i(r,r);
  double d, log_d;

  sym_inverse(VAR, VAR_i, log_d, d, 1e-5);
  score = beta.transpose() * ( VAR_i * beta );
}

#endif
