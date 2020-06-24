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

  // Rcout << "BFGS" << endl;
  uword p = x.n_cols;
  Rcout.precision(20);
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

    // Rcout << "Iter :" << iter << " Norm : " << norm(g) << endl;

    if(norm(g)<tolerance) break;

    mat step = -h*g;

    //Backtracking
    double t = wolfe_line_search(x,y,rho,u,z,step);

    // Rcout << " t : " << t << " Primal : " << primal(y,lin_pred)+0.5*rho*pow(norm(z-u),2) ;

    mat dz = t*step;
    mat dgrad = gradient(x,y,x*(z+dz))+rho*(z+dz-u) - g;
    
    // Rcout << " 1/ETA : " << dot(dg,dx) << " ETA :" << 1/dot(dg,dx)  << endl;
    double eta = 1.0/dot(dgrad,dz);

    

    if(dot(dgrad,dz)==0){
      Rcout << "Dot 0 " << endl;

      Rcout << std::fixed << g << endl;
      Rcout << "---------" << endl;
      Rcout << std::fixed << gradient(x,y,x*(z+dz))+rho*(z+dz-u) << endl;
      Rcout << "---------" << endl;
      dz.print();
      Rcout << "---------" << endl;
      dgrad.print();
      Rcout << "---------" << endl;
      Rcout << norm(g);
      // exit(0);
      break;
    }

    // Inverse Hessian update
    h = (I-eta*dz*dgrad.t())*h*(I-eta*dgrad*dz.t()) + eta*dz*dz.t();

    
    z = z + t*step;

    // if(i%20 == 0){
    //   Rcout << "z : ";
    //   z.print();
      // Rcout << "H inv : ";
      // h.print();
    // }

    Rcpp::checkUserInterrupt();
  }

  if(verbosity>=3){
    Rcout << " BFGS iterations : " << iter << " ";
  }

  // Rcout << endl << "BFGS ends" << endl;
  return z;
}
