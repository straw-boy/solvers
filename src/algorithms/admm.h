#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "../prox.h"
#include "../families/family.h"
#include "newton_raphson.h"

using namespace Rcpp;
using namespace arma;

// FISTA implementation
template <typename T>
Results Family::fitADMM(const T& x, const mat& y, vec lambda, double rho){
  uword p = x.n_cols;
  uword n = x.n_rows;
  uword m = y.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat z(p,m,fill::zeros);
  mat u(z);
  mat beta(z);
  mat beta_prev(z);
  
  // Ignore for now
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> time;


  uword passes=0;
  double alpha = 1.5;
  
  while(passes<max_passes){
    passes++;

    beta_prev = beta;
    
    beta = newton_raphson(x,y,rho,z-u,false);

    mat z_old = z;
    mat beta_hat = alpha*beta + (1 - alpha)*z_old;

    z = beta_hat + u;

    z.tail_rows(p_rows) = prox(z.tail_rows(p_rows), lambda/rho);

    u += (beta_hat-z);

    double r_norm = norm(beta - z);
    double s_norm = norm(rho*(z - z_old));

    double eps_primal = std::sqrt(n)*tol_abs + tol_rel*std::max(norm(beta), norm(z));
    double eps_dual = std::sqrt(n)*tol_abs + tol_rel*norm(rho*u);

    if (r_norm < eps_primal && s_norm < eps_dual)
        break;

    
    Rcpp::checkUserInterrupt();
  }

  double deviance = 2*primal(y, x*beta);

  Results res{beta,
            passes,
            primals,
            duals,
            time,
            deviance};

  return res;
}
