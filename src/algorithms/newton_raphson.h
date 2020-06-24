#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// // Newton-Raphson's method for optimization
// // Solves argmin_z { f(z) + rho/2*||z-u||^2 }
template <typename T>
mat Family::newton_raphson(const T& x, const mat&y, const double rho, const mat& u){

  Rcout.precision(4);
  uword p = x.n_cols;

  mat z(u);
  mat g(u);
  mat h(p,p);

  uword max_iter = 50;
  double tolerance = 1e-12;
  double alpha = 0.1;
  double gamma = 0.5;

  uword iter;

  for(iter = 0; iter < max_iter; iter++){
    
    mat lin_pred = x*z;
    
    g = gradient(x,y,lin_pred) + rho*(z-u);
    h = hessian(x,y,lin_pred);


    h.diag() += rho;

    mat step = -solve(h,g);
    double decrement = -dot(g,step);

    if( 0.5*decrement < tolerance){
      break; 
    }

    //Backtracking
    double t = 1.0;

    double f = primal(y,lin_pred)+(rho)/(2.0)*pow(norm(z-u),2);

    while( primal(y,x*(z+t*step))+(rho)/(2.0)*pow(norm(z+t*step-u),2)
            > (f + alpha*t*decrement) ){
        t = gamma*t;
        Rcpp::checkUserInterrupt();
    }

    z = z + t*step;

    Rcpp::checkUserInterrupt();
  }

  if(verbosity>=3){
    Rcout << "Newton-Raphson iterations " << iter << " ";
  }

  return z;
}
