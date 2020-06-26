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

    // Rcout << "Iter : " << iter << endl;
    // step.print();

    //Backtracking
    double t = wolfe_line_search(x,y,rho,u,z,step);
    // Rcout << "t is " << t << endl;

    mat dz = t*step;
    mat dgrad = gradient(x,y,x*(z+dz))+rho*(z+dz-u) - g;
    
    double eta = 1.0/dot(dgrad,dz);

    // // Does t actually satisfy wolfe conditions?
    // if( primal(y,x*(z+t*step))+0.5*rho*pow(norm(z+t*step-u),2) > primal(y,x*z)+0.5*rho*pow(norm(z-u),2) + 0.0001*t*dot(g,step)  
    //     || std::abs(dot(gradient(x,y,x*(z+t*step))+rho*(z+t*step-u),step)) > std::abs(0.9*dot(g,step)) ){
         
    //     Rcout << "Wolfe line search incorrect !!!!! " << t << endl;

    //     Rcout << "Dot: " << dot(dgrad,dz) << endl;

    //     // Rcout << std::fixed << g << endl;
    //     // Rcout << "---------" << endl;
    //     // Rcout << std::fixed << gradient(x,y,x*(z+dz))+rho*(z+dz-u) << endl;
    //     // Rcout << "---------" << endl;
    //     // dz.print();
    //     // Rcout << "---------" << endl;
    //     // dgrad.print();
    //     // Rcout << "---------" << endl;
    //     // Rcout << norm(g) << " " << t << " " <<  (-0.5)*dot(step,g) << endl;
    //     // break;
    // }
    // else if(dot(dgrad,dz)==0){
    //   Rcout << "DOT 0 but with wolfe" << endl;
    //   exit(0);
    // }

    // Change is too small so terminate
    if(dot(dgrad,dz)==0){
      break;
    }

    // Rcout << "dgrad: " << norm(dgrad) << " dz: " << norm(dz) << " eta: " << eta << endl;  

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
