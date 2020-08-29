#pragma once

#include <RcppArmadillo.h>
#include "../results.h"

using namespace Rcpp;
using namespace arma;

class Family {
protected:
  const bool intercept;
  const bool diagnostics;
  const uword max_passes;
  const double tol_rel_gap;
  const double tol_infeas;
  const double tol_abs;
  const double tol_rel;
  const double tol_coef;
  const uword verbosity;
  mat gaussian_hessian;


public:
  Family(const bool intercept,
         const bool diagnostics,
         const uword max_passes,
         const double tol_rel_gap,
         const double tol_infeas,
         const double tol_abs,
         const double tol_rel,
         const double tol_coef,
         const uword verbosity)
    : intercept(intercept),
      diagnostics(diagnostics),
      max_passes(max_passes),
      tol_rel_gap(tol_rel_gap),
      tol_infeas(tol_infeas),
      tol_abs(tol_abs),
      tol_rel(tol_rel),
      tol_coef(tol_coef),
      verbosity(verbosity) {}

  virtual double primal(const mat& y, const mat& lin_pred) = 0;

  virtual double dual(const mat& y, const mat& lin_pred) = 0;

  // this is not really the true gradient; it needs to multiplied by X^T
  virtual mat pseudoGradient(const mat& y, const mat& lin_pred) = 0;

  // this is not really the true hessian; it needs to left multiplied by X^T and right multiplies by X
  virtual mat pseudoHessian(const mat& y, const mat& lin_pred) = 0;

  template <typename T>
  mat gradient(T& x, const mat& y, const mat& lin_pred)
  {
    return x.t() * pseudoGradient(y, lin_pred);
  }

  template <typename T>
  mat hessian(const T& x, const mat& y, const mat& lin_pred)
  {
    if (name() == "multinomial"){
      mat tmp;
      tmp = x.t();

      vec lp_max = max(lin_pred, 1);
      vec lse = trunc_log(exp(-lp_max) 
                + sum(trunc_exp(lin_pred.each_col() - lp_max), 1)) + lp_max;
      
      mat prob_pred = trunc_exp(lin_pred.each_col() - lse);

      uword p = x.n_cols;
      uword m = y.n_cols;
      uword n = y.n_rows;

      mat H(p*m, p*m, fill::zeros);

      for (uword i = 0; i < m; i++) {
        for (uword j = 0; j < m; j++) {
          for (uword k = 0; k < n; k++) {
            vec v = tmp.col(k);
            H(span(i*p, (i+1)*p -1), span(j*p, (j+1)*p -1)) +=
              prob_pred(k, i)*((i == j) - prob_pred(k, j))*(v*v.t());
          }
        }
      }
      
      return H;

    } else if (name() == "gaussian") {
      if (gaussian_hessian.is_empty()) {
        gaussian_hessian = x.t();
        gaussian_hessian *= x;
      }
      return gaussian_hessian;
    }
    vec activation = pseudoHessian(y, lin_pred);
    mat xTx;
    xTx = x;
    xTx.each_col() %= activation;
    xTx = x.t()*xTx;
    return xTx;
  }


  template <typename T>
  mat newtonRaphson(const T& x, const mat& y, const double rho, const mat& u);

  template <typename T>
  mat bfgs(const T& x, const mat& y, const double rho, const mat& u);

  template <typename T>
  mat lbfgs(const T& x, const mat& y, const double rho, const mat& u);

  template <typename T>
  mat optimizeApproximation(const T& x, const mat& y, const double rho, const mat& u, const std::string opt_algo)
  {
    if (opt_algo == "bfgs")
        return bfgs(x, y, rho, u);
    else if (opt_algo == "nr")
        return newtonRaphson(x, y, rho, u);
    else
        return lbfgs(x, y, rho, u);
  }

  template <typename T>
  double zoom(const T& x,
              const mat& y,
              const double rho,
              const mat& u, 
              const mat& z,
              const mat& d,
              double t_low,
              double t_high);

  template <typename T>
  double wolfeLineSearch(const T& x,
                         const mat& y,
                         const double rho,
                         const mat& u,
                         const mat& z,
                         const mat& d);


  virtual rowvec fitNullModel(const mat& y, const uword n_classes) = 0;

  virtual std::string name() = 0;

  template <typename T>
  Results fitFISTA(const T& x, const mat& y, vec lambda);

  template <typename T>
  Results fitADMM(const T& x, const mat& y, vec lambda, const std::string opt_algo, const double rho = 1.0);

  template <typename T>
  mat scaled_prox(const T& x, const vec& y, const mat& beta, const mat& H, const vec& lambda);

  template <typename T>
  Results fitProximalNewton(const T& x, const mat& y, vec lambda);

  template <typename T>
  Results fitProximalQuasiNewton(const T& x, const mat& y, vec lambda);

  template <typename T>
  Results fitPN(const T& x, const mat& y, vec lambda, const std::string hessian_calc){
    if (hessian_calc == "exact")
      return fitProximalNewton(x, y, lambda);
    else 
      return fitProximalQuasiNewton(x, y, lambda);
  }
  
};
