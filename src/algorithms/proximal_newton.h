#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "scaled_prox.h"
#include "../scaled_prox.h"
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// Proximal Newton algorithm
template <typename T>
Results Family::fitProximalNewton(const T& x, const mat& y, vec lambda)
{
  uword n = y.n_rows;
  uword p = x.n_cols;
  uword m = y.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat beta(p, m, fill::zeros);
  mat beta_tilde(beta);

  mat lin_pred(n, m);
  mat grad;
  mat hess;

  // line search parameters
  const double alpha = 0.5;
  const double eta = 0.5;

  vec hess_correction(p, fill::ones);

  vec activation(p, fill::ones);

  // diagnostics
  wall_clock timer;
  std::vector<double> loss;
  std::vector<double> time;

  if (diagnostics) {
    loss.reserve(max_passes);
    time.reserve(max_passes);
    timer.tic();
  }

  // main loop
  uword passes = 0;
  while (passes < max_passes) {
    ++passes;

    lin_pred = x*beta;
  
    double f = primal(y, lin_pred);
    double g = dot(sort(abs(vectorise(beta.tail_rows(p_rows))),
                        "descending"), lambda);
    double obj = f + g;

    if (diagnostics) {
        loss.push_back(obj);
        time.push_back(timer.toc());
    }
    
    if (verbosity >= 3) {
        Rcout << "pass: "         << passes
              << ", objective: "  << obj 
              << endl;
    }

    grad = gradient(x, y, lin_pred);
    hess = hessian(x, y, lin_pred);

    if (name() != "multinomial") {
      activation = pseudoHessian(y, lin_pred);
    }

    // if (rcond(hess) < 1e-16) {
    //   Rcout << "corrected hessian" << endl;
    //   hess.diag() += hess_correction;
    // }

    beta_tilde = beta - solve(hess, grad);
    // beta_tilde = beta - solve(hess, grad, solve_opts::fast);

    beta_tilde = scaled_prox(x, activation, beta_tilde, hess, lambda);
    
    mat d = beta_tilde - beta;

    // Backtracking line search
    double t = 1.0;
    double dTgrad = dot(d, grad);
    double obj_new;
    while (true) {
      mat beta_new = beta + t*d;

      double f_new = primal(y, x*(beta_new));
      double g_new = dot(sort(abs(vectorise(beta_new.tail_rows(p_rows))),
                              "descending"), lambda);
      obj_new = f_new + g_new;

      if (obj + alpha*(t*dTgrad + g_new - g)  >= obj_new) {
        break;
      } else {
        t *= eta;
      }
      checkUserInterrupt();
    }

    beta += t*d;

    if (norm(t*d) < tol_coef)
      break;
    
    hess_correction = abs(t*d);
    
    if (passes % 10 == 0)
      checkUserInterrupt();
    
  }

  double deviance = 2*primal(y, x*beta);

  // res.diagnosticsLoss contains 'objective value'
  // at index 0
  Results res{beta,
              passes,
              {loss},
              time,
              deviance};

  return res;
}
