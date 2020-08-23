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
  mat g(size(z));
  mat h(p*m, p*m);
  mat step(size(z));

  uword max_iter = 50;
  double tolerance = 1e-15;
  double alpha = 0.1;
  double gamma = 0.5;

  // bool use_woodbury 
  //   = (name() != "multinomial" && n < p)? true: false;

  bool use_woodbury = false;

  mat xxT;
  vec activation;

  // vec hessian_corr(n, fill::ones);
  // double err = 1.0;

  if (use_woodbury) {
    xxT = x*x.t();
  }

  uword iter = 0;

  while (iter < max_iter){
    
    mat lin_pred = x*z;
    
    g = gradient(x, y, lin_pred) + rho*(z-u);
    
    if (use_woodbury) {
      // activation = pseudoHessian(y, lin_pred);

      // mat tmp = xxT;
      // tmp.diag() += rho/activation;
      // double epsilon = 1e-8;
      // if (rcond(tmp) < 1e-16) {
      //   Rcout << "hessian problem " << rcond(tmp) << endl;    
      //   hessian_corr.print();
      //   tmp.diag() += hessian_corr;
      //   Rcout << "-------" << rcond(tmp) << endl;
      // }

      // step = x.t()*solve(tmp, x*g) - g;
      // step /= rho;

    } else {

      h = hessian(x, y, lin_pred);
      h.diag() += rho;
      step = -reshape(solve(h, vectorise(g)),size(z));

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

    // hessian_corr = 0*hessian_corr;
    // for(int l=0;l<n;l++)
    //   hessian_corr(l) = std::abs(t*step(l));
    // hessian_corr += (dot(t*step, gradient(x, y, lin_pred) + rho*(z-u) - g));
    
    Rcpp::checkUserInterrupt();
    
    iter++;
  }

  if (verbosity >= 3) {
    Rcout << "Newton-Raphson iterations " << iter << " ";
  }

  return z;
}
