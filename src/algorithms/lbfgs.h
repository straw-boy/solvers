#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"
#include "wolfe_line_search.h"

using namespace Rcpp;
using namespace arma;

// // LBFGS Algorithm.
// // Solves argmin_z { f(z) + rho/2*||z-u||^2 }
template <typename T>
mat Family::lbfgs(const T& x, const mat& y, const double rho, const mat& u)
{

  mat z(u);
  mat g(u);
  
  const int max_iter = 2000;
  const double tolerance = 1e-5;
  const int l = 20;

  std::deque<mat> dz;
  std::deque<mat> dgrad;
  std::deque<double> eta;

  std::vector<double> alpha(l);

  // Scaling parameter
  double gamma = 1.0;

  mat lin_pred = x*z;
  g = gradient(x, y, lin_pred) + rho*(z-u);

  int iter = 0;

  while (iter < max_iter) {

    if (sqrt(accu(square(g))) < tolerance)
      break;

    // Two Loop recursion begins
    mat q = g;
    for (int j = std::min(iter, l)-1; j >= 0; j--) {
      alpha[j] = eta[j]*accu(dz[j] % q);
      q = q - alpha[j]*dgrad[j];
    }

    q = gamma*q;

    for (int j = 0; j < std::min(iter, l); j++) {
      q += dz[j]*(alpha[j]-eta[j]*accu(dgrad[j] % q));
    }
    // Two Loop recursion ends

    mat step = -q;

    // Getting rid of older changes in gradient and change in z
    if (iter > l) {
      dz.pop_front();
      dgrad.pop_front();
      eta.pop_front();
    }

    // Line Search for step length
    double t = wolfeLineSearch(x, y, rho, u, z, step);
    
    // Store change in z and gradient
    dz.push_back(t*step);
    z += dz.back();
    dgrad.push_back(gradient(x, y, x*z) + rho*(z-u) - g);

    double eta_inv = accu(dgrad.back() % dz.back());
    if (eta_inv == 0)
      break;
    eta.push_back(1.0/eta_inv);

    // Updating scaling parameter for next iteration
    gamma = eta_inv / accu(dgrad.back() % dgrad.back());

    g += dgrad.back();

    Rcpp::checkUserInterrupt();
    
    iter++;
  }

  if (verbosity >= 3) {
    Rcout << " L-BFGS iterations : " << iter << " ";
  }
  
  return z;
}
