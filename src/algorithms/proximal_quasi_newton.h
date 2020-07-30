#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "../lbfgs_utils.h"
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// Proximal Quasi-Newton algorithm. 
// Uses L-BFGS for Hessian approximation.
template <typename T>
Results Family::fitProximalQuasiNewton(const T& x, const mat& y, vec lambda)
{
  
  uword p = x.n_cols;
  uword m = y.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat beta(p, m, fill::zeros);
  mat beta_tilde(beta);

  // line search parameters
  const double alpha = 0.5;
  const double eta = 0.5;

  // diagnostics
  wall_clock timer;
  std::vector<double> loss;
  std::vector<double> time;

  if (diagnostics) {
    loss.reserve(max_passes);
    time.reserve(max_passes);
    timer.tic();
  }

  LBFGS lbfgs;
  
  mat lin_pred = x*beta;
  mat grad = gradient(x, y, lin_pred);

  // main loop
  uword passes = 0;
  while (passes < max_passes) {
    ++passes;
  
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

    beta_tilde = beta - lbfgs.inverseHessianProduct(grad);

    beta_tilde = lbfgs.scaled_prox(beta_tilde, lambda);
    
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
    
    beta_tilde = beta + t*d;

    lin_pred = x*beta_tilde;

    mat grad_new = gradient(x, y, lin_pred);

    if (lbfgs.updateParams(beta_tilde - beta, grad_new - grad)){
      break;
    }

    if (norm(t*d) < tol_coef)
      break;
    
    grad = grad_new;
    beta = beta_tilde;

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

