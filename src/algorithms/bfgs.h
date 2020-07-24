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
  uword m = y.n_cols;
  mat z(u);
  mat g(u);
  mat I(p*m, p*m, fill::eye);

  mat h(I);

  uword max_iter = 1000;
  double tolerance = 1e-6;

  uword iter = 0;

  while (iter < max_iter) {

    mat lin_pred = x*z;
    g = gradient(x, y, lin_pred) + rho*(z-u);

    if (sqrt(accu(square(g))) < tolerance)
      break;

    mat step = -reshape(h*vectorise(g.t()),size(z.t())).t();

    // Line Search for step length
    double t = wolfeLineSearch(x, y, rho, u, z, step);

    mat dz = t*step;
    mat dgrad = gradient(x, y, x*(z+dz)) + rho*(z+dz-u) - g;

    // Change is too small so terminate
    if (accu(dgrad % dz) == 0)
      break;

    double eta = 1.0/accu(dgrad % dz);

    // Inverse Hessian update
    vec dz_vec = vectorise(dz.t());
    vec dgrad_vec = vectorise(dgrad.t());
    h = (I-eta*dz_vec*dgrad_vec.t())*h*(I-eta*dgrad_vec*dz_vec.t()) + eta*dz_vec*dz_vec.t();
    
    z = z + t*step;

    Rcpp::checkUserInterrupt();

    iter++;
  }

  if (verbosity >= 3) {
    Rcout << " BFGS iterations : " << iter << " ";
  }

  return z;
}