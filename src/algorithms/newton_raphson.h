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
  uword n = x.n_rows;
  uword m = y.n_cols;

  mat z(u);
  mat g(u);
  mat h(p*m, p*m);

  uword max_iter = 50;
  double tolerance = 1e-12;
  double alpha = 0.1;
  double gamma = 0.5;

  uword iter = 0;

  while (iter < max_iter){
    
    mat lin_pred = x*z;
    
    g = gradient(x, y, lin_pred) + rho*(z-u);

    mat step;

    if (n >= p || name() == "multinomial") {
      h = hessian(x, y, lin_pred);
      h.diag() += rho;
      step = -reshape(solve(h,vectorise(g)),size(z));
    } else {
      // vec activation = pseudoHessian(y, lin_pred);
    
      // step = x*g;
      // step /= rho;

      // diagmat(activation).print();

      // mat tmp = diagmat(activation);
      // tmp = solve(tmp, x);
      // tmp = x*tmp;
      // tmp.diag() += 1/activation;

      // step = solve(tmp, step);
      // step = x.t()*step;
      // step.diag() -= 1;
      // step /= rho;
    }
    

    double decrement = -accu(g % step);
    
    if (0.5*decrement*decrement < tolerance) {
      break; 
    }

    //Backtracking
    double t = 1.0;

    double f = primal(y, lin_pred) + 0.5*rho*accu(square(z-u));

    while (primal(y, x*(z+t*step)) + 0.5*rho*accu(square(z+t*step-u))
            > (f + alpha*t*decrement)) {
        t = gamma*t;
        Rcpp::checkUserInterrupt();
    }

    z = z + t*step;
    
    Rcpp::checkUserInterrupt();
    
    iter++;
  }

  if (verbosity >= 3) {
    Rcout << "Newton-Raphson iterations " << iter << " ";
  }

  return z;
}
