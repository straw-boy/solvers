#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "../prox.h"
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
  mat grad(p, m, fill::zeros);
  mat hess(p, p, fill::zeros);

  // line search parameters
  const double alpha = 0.5;
  const double eta = 0.5;

  const double tol = 1e-10;

  // diagnostics
  wall_clock timer;
  std::vector<double> loss;
  std::vector<double> time;

  if(diagnostics){
    loss.reserve(max_passes);
    time.reserve(max_passes);
    timer.tic();
  }

  // main loop
  uword passes = 0;

  beta += 100;

  lambda.print();

  while (passes < max_passes) {
    ++passes;

    lin_pred = x*beta;

    double f = primal(y, lin_pred);
    double g = dot(sort(abs(vectorise(beta.tail_rows(p_rows))),
                        "descending"), lambda);
    double obj = f + g;

    grad = gradient(x, y, lin_pred);
    hess = hessian(x, y, lin_pred);

    if (verbosity >= 3) {
      Rcout << "pass: "         << passes
            << ", objective: "  << obj;
    }

    if (diagnostics) {
        loss.push_back(obj);
        time.push_back(timer.toc());
    }
    
    beta_tilde = beta - solve(hess, grad);

    beta_tilde = scaled_prox(beta_tilde, hess, lambda);

    mat d = beta_tilde - beta;

    if (verbosity >= 3) {
      Rcout << " , normd: " << norm(d);
    }

    // Backtracking line search
    double t = 1.0;
    double dTgrad = dot(d, grad);
    double f_new;
    double g_new;
    double obj_new;
    while (true) {
      mat beta_new = beta + t*d;
      f_new = primal(y, x*(beta_new));
      g_new = dot(sort(abs(vectorise(beta_new.tail_rows(p_rows))),
                         "descending"), lambda);
      obj_new = f_new + g_new;

      if (obj + alpha*(t*dTgrad + g_new - g)  >= obj_new) {
        break;
      } else {
        t *= eta;
      }
      checkUserInterrupt();
    }

    Rcout << ", t: " << t;
    
    beta += t*d;

    if (verbosity >= 3) {
      Rcout << " df: " << obj_new-obj << endl;
    }

    if (obj_new-obj < tol)
      break;

    if (passes % 100 == 0)
      checkUserInterrupt();
    
  }
  
  if (diagnostics) {
    loss.resize(passes);
    time.resize(passes);
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
