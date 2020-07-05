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
  uword n = x.n_rows;
  uword m = y.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat z(p, m, fill::zeros);
  mat u(z);
  mat beta(z);
  mat beta_prev(z);

  std::vector<double> loss;
  std::vector<double> eps_primals;
  std::vector<double> eps_duals;
  std::vector<double> time;

  wall_clock timer;

  if (diagnostics) {
    loss.reserve(max_passes);
    eps_primals.reserve(max_passes);
    eps_duals.reserve(max_passes);
    time.reserve(max_passes);
    timer.tic();
  }

  double alpha = 1.5;
  uword passes = 0;

  while (passes < max_passes) {
    passes++;

    beta = optimizeApproximation(x, y, rho, z-u, opt_algo);

    mat z_old = z;
    mat beta_hat = alpha*beta + (1 - alpha)*z_old;

    z = beta_hat + u;

    z.tail_rows(p_rows) = prox(z.tail_rows(p_rows), lambda/rho);

    u += (beta_hat-z);

    double r_norm = norm(beta - z);
    double s_norm = norm(rho*(z - z_old));

    double eps_primal = std::sqrt(n)*tol_abs + tol_rel*std::max(norm(beta), norm(z));
    double eps_dual = std::sqrt(n)*tol_abs + tol_rel*norm(rho*u);

    if (diagnostics) {
      loss.push_back(primal(y, x*beta) + 
                     dot(sort(abs(vectorise(beta.tail_rows(p_rows))),"descending"), lambda));
      eps_primals.push_back(r_norm);
      eps_duals.push_back(s_norm);
      time.push_back(timer.toc());
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

  if (diagnostics) {
    loss.resize(passes);
    eps_primals.resize(passes);
    eps_duals.resize(passes);
    time.resize(passes);
  }

  double deviance = 2*primal(y, x*beta);
  
  // res.diagnosticsLoss contains 'objective value', 'residual primals' 
  // and 'residual duals' at indices 0, 1 and 2 respectively
  Results res{beta,
              passes,
              {loss, eps_primals, eps_duals},
              time,
              deviance};

  return res;
}
