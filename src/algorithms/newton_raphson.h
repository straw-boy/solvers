#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// Newton-Raphson's method for optimization
// Solves argmin_z { f(z) + rho/2*||z-u||^2 }
template <typename T>
mat Family::newton_raphson(const T& x, const mat&y, const double rho, const mat& u, bool quasi){
  // Rcout << "Newton raphson's method starts" << endl;
  
  uword p = u.n_cols;

  mat z(u);
  mat g(u);
  mat h(p,p);

  uword max_iter = 500;
  double tolerance = 1e-8;


  for(uword i=0;i<max_iter;i++){
    
    mat lin_pred = x*z;
    
    // z.t().print();
    // Rcout << primal(y,lin_pred)+(rho)/(2.0)*norm(z-u,2) << endl;
    
    g = gradient(x,y,lin_pred)+rho*(z-u);
    h = hessian(x,y,lin_pred);

    h.diag() += rho;

    mat step = solve(h,g);
    // mat decrement = g.t()*step;

    // if( decrement < tolerance)
    //   break;

    // TODO - Backtracking

    z = z - step;
    if(step.is_zero(tolerance)){
      break;
    }
    
    Rcpp::checkUserInterrupt();
  }

  // mat ans = hessian(x,y,x*u);
  // ans.diag() += rho;
  // vec anns = solve(ans,x.t()*y+rho*u);
  // anns.print();
  // Rcout << primal(y,x*anns)+(rho)/(2.0)*norm(anns-u,2) << endl;

  // Rcout << "Newton raphson's method ends" << endl;
  return z;
}