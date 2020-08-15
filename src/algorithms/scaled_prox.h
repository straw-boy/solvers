#pragma once

#include <RcppArmadillo.h>
#include "../prox.h"

using namespace Rcpp;
using namespace arma;


// Solvers argmin_v { 0.5*(v-beta).t()*H*(v-beta) + J(v, lambda) }
// where J(v, lambda) is the slope penalty.
// This subproblem is solved using ADMM.
template <typename T>
mat Family::scaled_prox(const T& x, const vec& activation, const mat& beta, const mat& H, const vec& lambda)
{
  uword n = x.n_rows;
  uword p = beta.n_rows;
  uword m = beta.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat v(beta);
  mat z(v);
  mat u(v);

  mat xx;

  const double alpha = 1.5;
  const double rho = 1.0;
  const double tol_abs = 1e-6;
  const double tol_rel = 1e-5;

  mat H_beta = H*beta;

  if (n != p || name() == "multinomial") {
    xx = H;
    xx.diag() += rho;
  } else {
    xx = x;
    xx = xx * x.t();
    xx.diag() += rho/activation;
  }
  mat U = chol(xx);
  mat L = U.t();

  uword iter = 0;
  uword max_iter = 500000;

  while (iter < max_iter) {
    iter++;
    
    v = H_beta + rho*(z - u);
    if (true || name() == "multinomial") {
      v = solve(trimatu(U), solve(trimatl(L), v));
    } else {
      v = v - x.t() * solve(trimatu(U), solve(trimatl(L), x*v));
      v /= rho;
    }

    mat z_old = z;
    mat v_hat = alpha*v + (1 - alpha)*z_old;

    z = v_hat + u;

    z.tail_rows(p_rows) = prox(z.tail_rows(p_rows), lambda/rho);

    u += (v_hat - z);

    double r_norm = norm(v - z);
    double s_norm = norm(rho*(z - z_old));

    double eps_primal = std::sqrt(p)*tol_abs + tol_rel*std::max(norm(v), norm(z));
    double eps_dual = std::sqrt(p)*tol_abs + tol_rel*norm(rho*u);

    if (r_norm < eps_primal && s_norm < eps_dual)
        break;

    if (iter % 1000 == 0)
      Rcpp::checkUserInterrupt();
  }
  Rcout << "SP: " << iter << endl;
  return v;
}
