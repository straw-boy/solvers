#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "../prox.h"
#include "../scaled_prox.h"
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// Proximal Newton
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

  mat I(p, p, fill::eye);

  // line search parameters
  double alpha = 0.5;
  double eta = 0.5;

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

  while (passes < max_passes) {
    
    lin_pred = x*beta;

    double f = primal(y, lin_pred);
    double g = dot(sort(abs(vectorise(beta.tail_rows(p_rows))),
                        "descending"), lambda);

    grad = gradient(x, y, lin_pred);
    hess = hessian(x, y, lin_pred);
    
    if (verbosity >= 3) {
      Rcout << "pass: "         << passes
            << ", objective: "  << f + g
            << std::endl;
    }

    if (diagnostics) {
        loss.push_back(f+g);
        time.push_back(timer.toc());
        timer.tic();
    }
    
    beta_tilde = beta - solve(hess, grad);

    beta_tilde = scaled_prox(beta_tilde, hess, lambda);

    ///experiment
    mat a = scaled_prox(beta_tilde, I, lambda);
    mat b = prox(beta_tilde, lambda);

    if (!(a-b).is_zero(1e-3)) {
      Rcout << "Messed up scaled prox" << endl;
      (a-b).print();
    }
    ///

    mat d = beta_tilde - beta;
    
    // Backtracking line search
    double t = 1.0;
    double dTgrad = dot(d,grad);
    while (true) {
      vec beta_new = beta + t*d;
      double f_new = primal(y, x*(beta_new));
      double g_new = dot(sort(abs(vectorise(beta_new.tail_rows(p_rows))),
                              "descending"), lambda);

      if (f + alpha*(t*dTgrad + g_new - g)  >= f_new*(1 - 1e-12)) {
        break;
      } else {
        t *= eta;
      }
      checkUserInterrupt();
    }
    Rcout << " t : " << t << endl;

    beta += t*d;

    ++passes;
    
    double df = primal(y, x*beta) + 
                dot(sort(abs(vectorise(beta.tail_rows(p_rows))),
                         "descending"), lambda);

    df = df - f - g;
    // Rcout << "df: " << df << endl;
    if (std::abs(df) < 1e-10)
      break; 

    
    if (passes % 100 == 0)
      checkUserInterrupt();
    
  }

  Rcout << "Final obj: " << primal(y, x*beta) + dot(sort(abs(vectorise(beta.tail_rows(p_rows))),
                                                          "descending"), lambda) << endl;

  loss.reserve(max_passes);
  time.reserve(max_passes);

  double deviance = 2*primal(y, lin_pred);

  Results res{beta,
              passes,
              loss,
              loss,
              time,
              deviance};

  return res;
}
