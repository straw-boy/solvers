#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "../prox.h"
#include "../infeasibility.h"
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// FISTA implementation
template <typename T>
Results Family::fitFISTA(const T& x, const mat& y, vec lambda){
  
  uword n = y.n_rows;
  uword p = x.n_cols;
  uword m = y.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat beta(p, m, fill::zeros);
  beta += 0.1;
  mat beta_tilde(beta);
  mat beta_tilde_old(beta);

  mat lin_pred(n, m);
  mat grad(p, m, fill::zeros);

  double learning_rate = 1.0;

  // line search parameters
  double eta = 0.5;

  // FISTA parameters
  double t = 1;

  // diagnostics
  wall_clock timer;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> time;

  if(diagnostics){
    primals.reserve(max_passes);
    duals.reserve(max_passes);
    time.reserve(max_passes);
    timer.tic();
  }
  

  // main loop
  uword passes = 0;
  while (true) {
    
    lin_pred = x*beta;

    double g = primal(y, lin_pred);
    double h = dot(sort(abs(vectorise(beta.tail_rows(p_rows))),
                        "descending"), lambda);
    double f = g + h;
    double G = dual(y, lin_pred);

    grad = gradient(x, y, lin_pred);
    double infeas = lambda.n_elem > 0 ? infeasibility(grad.tail_rows(p_rows), lambda) : 0.0;

    double small = std::sqrt(datum::eps);

    bool optimal =
      (std::abs(f - G)/std::max(small, std::abs(f)) < tol_rel_gap);

    bool feasible = 
      lambda.n_elem > 0 ? infeas <= std::max(small, tol_infeas*lambda(0)) : true;

    
    if (verbosity >= 3) {
      Rcout << "objective " << f << " ";
      Rcout << "pass: "            << passes
            << ", duality-gap: "   << std::abs(f - G)/std::abs(f)
            << ", infeasibility: " << infeas
            << std::endl;
    }

    if (diagnostics) {
        primals.push_back(f);
        duals.push_back(G);
        time.push_back(timer.toc());
        timer.tic();
    }

    // Rcout << "--------" << endl;
    // beta.print();
    // Rcout << "--------" << endl;
    // grad.print();
    // Rcout << "--------" << endl;


    if (optimal && feasible){
      break;
    }

    beta_tilde_old = beta_tilde;

    double g_old = g;
    double t_old = t;

    // Backtracking line search
    while (true) {
      // Update coefficients
      beta_tilde = beta - learning_rate*grad;

      beta_tilde.tail_rows(p_rows) =
        prox(beta_tilde.tail_rows(p_rows), lambda*learning_rate);

      vec d = vectorise(beta_tilde - beta);

      lin_pred = x*beta_tilde;

      g = primal(y, lin_pred);

      double q = g_old
        + dot(d, vectorise(grad))
        + (1.0/(2*learning_rate))*accu(square(d));

        if (q >= g*(1 - 1e-12)) {
          break;
        } else {
          learning_rate *= eta;
        }
        checkUserInterrupt();
    }

    // FISTA step
    t = 0.5*(1.0 + std::sqrt(1.0 + 4.0*t_old*t_old));
    beta = beta_tilde + (t_old - 1.0)/t * (beta_tilde - beta_tilde_old);

    if (passes % 100 == 0)
      checkUserInterrupt();

    ++passes;
    
  }

  double deviance = 2*primal(y, lin_pred);

  Results res{beta,
              passes,
              primals,
              duals,
              time,
              deviance};

  return res;
}
