#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"
#include "wolfe_line_search.h"

using namespace Rcpp;
using namespace arma;

// // BFGS Algorithm.
// // Solves argmin_z { f(z) + (rho/2)*||z-u||^2 }
template <typename T>
mat Family::bfgs(const T& x, const mat& y, const double rho, const mat& u)
{
  
  uword p = x.n_cols;
  mat z(u);
  mat g(u);
  mat I(p, p, fill::eye);

  mat h(I);

  uword max_iter = 1000;
  double tolerance = 1e-6;

  uword iter = 0;

  while (iter < max_iter) {

    mat lin_pred = x*z;
    g = gradient(x, y, lin_pred) + rho*(z-u);

    if (norm(g) < tolerance)
      break;

    mat step = -h*g;

    // Line Search for step length
    double t = wolfeLineSearch(x, y, rho, u, z, step);

    mat dz = t*step;
    mat dgrad = gradient(x, y, x*(z+dz)) + rho*(z+dz-u) - g;

    // Change is too small so terminate
    if (dot(dgrad, dz) == 0)
      break;

    double eta = 1.0/dot(dgrad, dz);

    // Inverse Hessian update
    h = (I-eta*dz*dgrad.t())*h*(I-eta*dgrad*dz.t()) + eta*dz*dz.t();
    
    z = z + t*step;

    Rcpp::checkUserInterrupt();

    iter++;
  }

  if (verbosity >= 3) {
    Rcout << " BFGS iterations : " << iter << " ";
  }

  return z;
}
