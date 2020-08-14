#pragma once

#include <RcppArmadillo.h>
#include "prox.h"

using namespace Rcpp;
using namespace arma;


// Solvers argmin_x { 0.5*(x-beta).t()*H*(x-beta) + J(x,lambda) }
// where J(x,lambda) is the slope penalty.
// This subproblem is solved using ADMM.
inline mat scaled_prox2(const mat& beta, const mat& H, const vec& lambda)
{
  uword p = beta.n_rows;
  uword m = beta.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat x(beta);
  mat z(x);
  mat u(x);

  mat I(p, p, fill::eye);

  const double alpha = 1.5;
  const double rho = 1.0;
  const double tol_abs = 1e-6;
  const double tol_rel = 1e-5;

  uword iter = 0;
  uword max_iter = 500000;

  while (iter < max_iter) {
    iter++;
    
    x = solve(H + rho*I, H*beta + rho*(z - u));
    
    mat z_old = z;
    mat x_hat = alpha*x + (1 - alpha)*z_old;

    z = x_hat + u;

    z.tail_rows(p_rows) = prox(z.tail_rows(p_rows), lambda/rho);

    u += (x_hat-z);

    double r_norm = norm(x - z);
    double s_norm = norm(rho*(z - z_old));

    double eps_primal = std::sqrt(p)*tol_abs + tol_rel*std::max(norm(x), norm(z));
    double eps_dual = std::sqrt(p)*tol_abs + tol_rel*norm(rho*u);

    if (r_norm < eps_primal && s_norm < eps_dual)
        break;

    if (iter % 1000 == 0)
      Rcpp::checkUserInterrupt();
  }
  
  return x;
}
