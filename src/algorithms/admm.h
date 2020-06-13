#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "../utils.h"
#include "../prox.h"
#include "../families/family.h"
#include "newton_raphson.h"

using namespace Rcpp;
using namespace arma;

// FISTA implementation
template <typename T>
Results Family::fitADMM(const T& x, const mat& y, vec lambda, double rho){
  Rcout << "ADMM starts" << endl;

  uword p = x.n_cols;
  uword m = y.n_cols;

  mat z(p,m,fill::zeros);
  mat u(z);
  mat beta(z);
  mat beta_prev(z);

  // Ignore for now
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> time;


  uword passes=0;

  while(passes<max_passes){
    passes++;

    beta_prev = beta;
    beta = newton_raphson(x,y,rho,z-u,false);

    z = beta+u;
    z.tail_rows(lambda.n_elem) = prox(z.tail_rows(lambda.n_elem), lambda/rho);

    u += (beta-z);

    if( (beta-beta_prev).is_zero()){
      Rcout << "ADMM Passes = " << passes << endl;
      break;
    }

  }

  double deviance = 2*primal(y, x*beta);
  // double deviance = 0.0;

  Results res{beta,
            passes,
            primals,
            duals,
            time,
            deviance};

  return res;
}