#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// // Newton-Raphson's method for optimization
// // Solves argmin_z { f(z) + rho/2*||z-u||^2 }
template <typename T>
mat Family::newton_raphson(const T& x, const mat&y, const double rho, const mat& u, bool quasi){

  uword p = x.n_cols;

  mat z(u);
  mat g(u);
  mat h(p,p);

  uword max_iter = 20;
  double tolerance = 1e-8;
  double alpha = 0.2;
  double gamma = 0.5;

  uword i;

  for(i=0;i<max_iter;i++){
    
    mat lin_pred = x*z;
    
    g = gradient(x,y,lin_pred)+rho*(z-u);
    h = hessian(x,y,lin_pred);

    h.diag() += rho;

    mat step = solve(h,g);
    double decrement = dot(g,step);

    // Rcout << "Decrement:" << decrement << "  Newton Passes: " << i << endl;
    // g.t().print();
    // step.t().print();
    
    if( decrement < tolerance){
      break; 
    }

    //Backtracking; Doesn't work right now.
    double t = 1.0;

    // double f = primal(y,lin_pred)+(rho)/(2.0)*norm(z-u,2);
    
    // while( primal(y,x*(z-t*step))+(rho)/(2.0)*norm(z-t*step-u,2) 
    //         > (f - alpha*t*decrement) ){
    //     t = gamma*t;
    //     Rcpp::checkUserInterrupt();
    // }
            
    z = z - t*step;

    Rcpp::checkUserInterrupt();
  }

  mat xTx ;
  xTx = x.t()*x;
  xTx.diag() += rho;

  // if( !(z-solve(xTx,x.t()*y+rho*(u))).is_zero(1e-4) ){
  //   Rcout << "Messed up " << endl;
  //   z.print();
  //   Rcout << "----------" << endl;
  //   solve(xTx,x.t()*y+rho*(z-u)).print();
  //   Rcout << "----------" << endl;
  //   (z-solve(xTx,x.t()*y+rho*(u))).print();
  //   Rcout << "----------" << endl;
  //   u.print();
  //   Rcout << "----------" << endl;
  // }

  // Rcout << " Total Newton Passes: " << i << endl;

  return z;
}