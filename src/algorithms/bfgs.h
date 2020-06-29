#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"
#include "wolfe_line_search.h"

using namespace Rcpp;
using namespace arma;

// // BFGS Algorithm.
// // Solves argmin_z { f(z) + rho/2*||z-u||^2 }
template <typename T>
mat Family::bfgs(const T& x, const mat&y, const double rho, const mat& u){
  uword p = x.n_cols;
  mat z(u);
  mat g(u);
  mat I(p,p,fill::eye);

  mat h(I);

  uword max_iter = 1000;
  double tolerance = 1e-6;

  uword iter;

  for(iter = 0; iter < max_iter; iter++){

    mat lin_pred = x*z;
    g = gradient(x,y,lin_pred) + rho*(z-u);

    if(norm(g)<tolerance) break;

    mat step = -h*g;

    //Backtracking
    double t = wolfe_line_search(x,y,rho,u,z,step);
    // Rcout << "t is " << t << endl;

    mat dz = t*step;
    mat dgrad = gradient(x,y,x*(z+dz))+rho*(z+dz-u) - g;
    
    double eta = 1.0/dot(dgrad,dz);

    // Change is too small so terminate
    if(dot(dgrad,dz)==0){
      break;
    }

    // Inverse Hessian update
    h = (I-eta*dz*dgrad.t())*h*(I-eta*dgrad*dz.t()) + eta*dz*dz.t();
    
    z = z + t*step;

    Rcpp::checkUserInterrupt();
  }

  if(verbosity>=3){
    Rcout << " BFGS iterations : " << iter << " ";
  }

  return z;
}
