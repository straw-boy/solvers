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
  const uword verbosity;

public:
  Family(const bool intercept,
         const bool diagnostics,
         const uword max_passes,
         const double tol_rel_gap,
         const double tol_infeas,
         const double tol_abs,
         const double tol_rel,
         const uword verbosity)
    : intercept(intercept),
      diagnostics(diagnostics),
      max_passes(max_passes),
      tol_rel_gap(tol_rel_gap),
      tol_infeas(tol_infeas),
      tol_abs(tol_abs),
      tol_rel(tol_rel),
      verbosity(verbosity) {}

  virtual double primal(const mat& y, const mat& lin_pred) = 0;

  virtual double dual(const mat& y, const mat& lin_pred) = 0;

  // this is not really the true gradient; it needs to multiplied by X^T
  virtual mat pseudoGradient(const mat& y, const mat& lin_pred) = 0;

  // this is not really the true hessian; it needs to left multiplied by X^T and right multiplies by X
  virtual mat pseudoHessian(const mat& y, const mat& lin_pred) = 0;

  template <typename T>
  mat gradient(T& x, const mat&y, const mat& lin_pred)
  {
    return x.t() * pseudoGradient(y, lin_pred);
  }

  template <typename T>
  mat hessian(const T& x, const mat&y, const mat& lin_pred)
  {
    vec activation = pseudoHessian(y, lin_pred);
    mat xTx;
    xTx = x.t() * diagmat(activation) * x;
    return xTx;
  }

  template <typename T>
  mat newton_raphson(const T& x, const mat&y, const double rho, const mat& u, bool quasi = true);

  template <typename T>
  mat newton_raphsonq(const T& x, const mat&y, const double rho, const mat& u, bool quasi = true);

  virtual rowvec fitNullModel(const mat& y, const uword n_classes) = 0;

  virtual std::string name() = 0;

  template <typename T>
  Results fitFISTA(const T& x, const mat& y, vec lambda);

  template <typename T>
  Results fitADMM(const T& x, const mat& y, vec lambda, double rho = 1.0);
  
};
