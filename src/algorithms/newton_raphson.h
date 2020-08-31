#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// // Newton-Raphson's method for optimization
// // Solves argmin_z { f(z) + (rho/2)*||z-u||^2 }
template <typename T>
mat Family::newtonRaphson(const T& x, const mat& y, const double rho, const mat& u)
{

  Rcout.precision(4);
  uword p = x.n_cols;
  uword m = y.n_cols;

  mat z(u);
  mat g(size(z));
  mat h(p*m, p*m);
  mat step(size(z));

  uword max_iter = 50;
  double tolerance = 1e-15;

  // bool use_woodbury 
  //   = (name() != "multinomial" && n < p)? true: false;

  bool use_woodbury = false;

  mat xxT;
  vec activation;
  
  if (use_woodbury) {
    xxT = x*x.t();
  }

  uword iter = 0;

  while (iter < max_iter){
    
    mat lin_pred = x*z;
    
    g = gradient(x, y, lin_pred) + rho*(z-u);
    
    if (use_woodbury) {
      // This section is never executed as use_woodbury is false
      activation = pseudoHessian(y, lin_pred);

      mat tmp = xxT;
      tmp.diag() += rho/activation;
      step = x.t()*solve(tmp, x*g) - g;
      step /= rho;

    } else {

      h = hessian(x, y, lin_pred);
      h.diag() += rho;
      step = -reshape(solve(h, vectorise(g)),size(z));

    }

    double decrement = -accu(g % step);
    
    if (0.5*decrement*decrement < tolerance) {
      break; 
    }

    // Line Search for step length
    double t = wolfeLineSearch(x, y, rho, u, z, step);

    z = z + t*step;
    
    Rcpp::checkUserInterrupt();
    
    iter++;
  }

  if (verbosity >= 3) {
    Rcout << "Newton-Raphson iterations " << iter << " ";
  }

  return z;
}
