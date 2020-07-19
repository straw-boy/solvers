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
  
  uword n = y.n_rows;
  uword p = x.n_cols;
  uword m = y.n_cols;
  uword pmi = lambda.n_elem;
  uword p_rows = pmi/m;

  mat beta(p, m, fill::zeros);
  mat beta_tilde(beta);

  mat lin_pred(n, m);
  mat grad(p, m, fill::zeros);

  // mat hess(p, p, fill::zeros);

  // hess.diag() += 1;

  // mat I = hess;

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

  LBFGS lbfgs(lambda);
  // beta(0) = -2.398615e+00;
  lin_pred = x*beta;
  grad = gradient(x, y, lin_pred);

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
    // Rcout << "H check" << endl;
    // (hess-lbfgs.computeH(p)).print();
    // Rcout << "HV check" << endl;
    // (lbfgs.computeHv(grad)- hess*grad).print();
    // Rcout << "----------" << endl;
    // (lbfgs.computeHv2(grad)-hess*grad).print();
    // Rcout << "BV check" << endl;
    // (lbfgs.computeBv(grad)-inv(hess)*grad).print();
    // Rcout << "----------" << endl;
    // (lbfgs.computeBv2(grad)-solve(hess,grad)).print();
    // Rcout << "----------" << endl;
    

    beta_tilde = beta - lbfgs.computeHv(grad);

    beta_tilde = lbfgs.scaled_prox(beta_tilde);
    
    mat d = beta_tilde - beta;

    Rcout << "d is " << d(0) << endl;

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
    Rcout << " QPN t is : " << t << " norm(t*d): " << norm(t*d) << endl;
    beta_tilde = beta + t*d;
    // if (intercept)
    //   beta_tilde(0) = beta_tilde(0) - t*d(0) + d(0);
    
    lin_pred = x*beta_tilde;

    mat grad_new = gradient(x, y, lin_pred);

    if (lbfgs.updateParams(beta_tilde - beta, grad_new - grad)){
      Rcout << "UPDATE ENDED" << endl;
      // exit(0);
      break;
    }

    if (norm(t*d) < tol)
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
