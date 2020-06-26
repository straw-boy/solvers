#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"
#include "wolfe_line_search.h"

using namespace Rcpp;
using namespace arma;

// // LBFGS Algorithm.
// // Solves argmin_z { f(z) + rho/2*||z-u||^2 }
template <typename T>
mat Family::lbfgs(const T& x, const mat&y, const double rho, const mat& u){

  // Rcout << "LBFGS starts " << endl;
  uword p = x.n_cols;
  mat z(u);
  mat g(u);
  
  const int max_iter = 1000;
  const double tolerance = 1e-6;
  const int l = 50;

  std::deque <mat> dz;
  std::deque <mat> dgrad;
  std::deque <double> eta;

  std::vector <double> alpha(l);

  double gamma = 1.0;

  mat lin_pred = x*z;
  g = gradient(x,y,lin_pred) + rho*(z-u);

  int iter;

  for(iter = 0; iter < max_iter; iter++){

    if(norm(g) < tolerance) break;

    // Two Loop recursion
    mat q = g;
    for(int j = std::min(iter,l)-1; j >= 0; j--){
      alpha[j] = eta[j]*dot(dz[j],q);
      q = q - alpha[j]*dgrad[j];
    }

    q = gamma*q;

    for(int j = 0; j < std::min(iter,l); j++){
      q += dz[j]*(alpha[j]-eta[j]*dot(dgrad[j],q));
    }

    mat step = -q;
    
    if(iter > l){
      dz.pop_front();
      dgrad.pop_front();
      eta.pop_front();
    }

    // Backtracking 
    double t = wolfe_line_search(x,y,rho,u,z,step);
    Rcout << "Iter: " << iter << " t is " << t << " ";
    
    dz.push_back(t*step);
    z += dz.back();

    dgrad.push_back(gradient(x,y,x*z) + rho*(z-u) - g);

    double eta_inv = dot(dgrad.back(),dz.back());
    if(eta_inv == 0)
      break;
    eta.push_back(1.0/eta_inv);

    gamma = eta_inv / norm(dgrad.back());
    g += dgrad.back();

    Rcout << "norm : " << norm(g) << " dgrad: " << norm(dgrad.back()) << " dz: " << norm(dz.back()) << " eta: " << eta.back() << " eta_inv: " << eta_inv << endl;  

    Rcpp::checkUserInterrupt();
  }

  if(verbosity>=3){
    Rcout << " L-BFGS iterations : " << iter << " ";
  }
  // Rcout << "LBFGS ends " << endl;
  return z;
}
