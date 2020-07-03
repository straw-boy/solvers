#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "../prox.h"
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;


// ADMM implementation
template <typename T>
Results Family::fitADMM(const T& x, const mat& y, vec lambda, const std::string opt_algo, const double rho)
{
  
  uword p = x.n_cols;
  uword m = y.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat z(p, m, fill::zeros);
  mat u(z);
  mat beta(z);

  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> time;

  wall_clock timer;
  if (diagnostics)
      timer.tic();

  double alpha = 1.5;
  uword passes = 0;

  while (passes < max_passes) {
    passes++;

    Rcout << " obj: " << primal(y, x*beta) + dot(sort(abs(vectorise(beta.tail_rows(p_rows))),
                                  "descending"), lambda) << endl; 

    beta = optimizeApproximation(x, y, rho, z-u, opt_algo);

    mat z_old = z;
    mat beta_hat = alpha*beta + (1 - alpha)*z_old;

    z = beta_hat + u;

    z.tail_rows(p_rows) = prox(z.tail_rows(p_rows), lambda/rho);

    u += (beta_hat-z);

    double r_norm = norm(beta - z);
    double s_norm = norm(rho*(z - z_old));

    double eps_primal = std::sqrt(p)*tol_abs + tol_rel*std::max(norm(beta), norm(z));
    double eps_dual = std::sqrt(p)*tol_abs + tol_rel*norm(rho*u);

    if (diagnostics) {
      primals.push_back(r_norm);
      duals.push_back(s_norm);
      time.push_back(timer.toc());
      timer.tic();
    }

    if (verbosity >= 3) {
      Rcout << "pass: "              << passes
            << ", primal residual: " << r_norm
            << ", dual residual: "   << s_norm
            << std::endl;
    }

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
