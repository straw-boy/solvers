#pragma once

#include <RcppArmadillo.h>
#include "prox.h"
#include "infeasibility.h"

using namespace Rcpp;
using namespace arma;


class LBFGS {
private:
  // Number of variable and gradient changes
  // to be used
  const uword m = 7;
  // Scaling parameter of hessian
  double sigma = 1.0;
  // Scaling parameter of inverse hessian
  double gamma = 1.0;
  
  // Matrix storing changes in variable
  mat S;
  // Matrix storing changes in gradient
  mat Y;
  // Strictly lower traingle of S.t()*Y
  mat L;
  // Upper traingle of S.t()*Y
  mat R;
  // Diagonal of S.t()*Y
  mat D;

  // Utility matrices used in hessianProduct
  mat Q_Bv;
  mat A_Bv;

  // Utility matrices used in inverseHessianProduct
  mat Q_Hv;
  mat A_Hv;

  // Scaled Proximal parameters
  const double tol_rel_gap = 1e-5;
  const double tol_infeas = 1e-3;
  const double tol_rel_coef_change = 1e-3;
  uword max_passes = 10000;

public:
  // Updates S,Y,L,R,D.
  // Returns true to terminate main PQN loop, false otherwise.
  bool updateParams(const mat& s, const mat& y) {

    if (accu(s % y) == 0)
      return true;
    
    S.insert_cols(S.n_cols, vectorise(s));
    Y.insert_cols(Y.n_cols, vectorise(y));

    if (S.n_cols > m) {
      S.shed_col(0);
      Y.shed_col(0);
    }

    gamma = accu(s % y) / accu(y % y);
    sigma = 1/gamma;
    
    mat sTy = S.t()*Y;
    
    R = trimatu(sTy);
    L = sTy - R;
    D = diagmat(sTy);

    Q_Hv = join_rows(S, gamma*Y);
    Q_Bv = join_rows(sigma*S, Y);

    mat tmp;

    mat Rinv = inv(R);
    tmp = Rinv.t()*(D + gamma*Y.t()*Y)*Rinv;
    tmp = join_rows(tmp, -Rinv.t());
    A_Hv = join_cols(tmp, join_rows(-Rinv, zeros(size(R))));

    tmp = sigma*S.t()*S;
    tmp = join_rows(tmp, L);
    A_Bv = join_cols(tmp, join_rows(L.t(), -D));

    return false;

  }

  // Computes H*v, where H is inverse hessian approximation
  mat inverseHessianProduct(const mat& v) {
    vec col_v = vectorise(v);
    if (S.n_cols == 0) {
      return gamma*col_v;
    }

    mat ans = Q_Hv.t()*col_v;
    ans = A_Hv*ans;
    ans = Q_Hv*ans;
    ans = gamma*col_v + ans;

    return ans;
  }


  // Computes B*v, where B is hessian approximation
  mat hessianProduct(const mat& v) {
    vec col_v = vectorise(v);
    if (S.n_cols == 0) {
      return sigma*col_v;
    }

    mat ans = Q_Bv.t()*col_v;
    ans = inv(A_Bv)*ans;
    ans = Q_Bv*ans;
    ans = sigma*col_v - ans;

    return ans;
  }

  // Solvers argmin_x { 0.5*(x-beta).t()*B*(x-beta) + J(x,lambda) }
  // where J(x,lambda) is the slope penalty and B is hessian approximation.
  // This subproblem is solved using FISTA.
  mat scaled_prox(const mat& beta, vec lambda) {
    uword p = beta.n_rows;
    uword m = beta.n_cols;
    uword pmi = lambda.n_elem;
    uword p_rows = pmi/m;

    mat x(beta);
    mat x_prev(x);
    mat x_tilde(x);
    mat x_tilde_old(x);

    double learning_rate = 1.0;
    
    double alpha = 0.5;

    double t = 1;

    uword passes = 0;
    while (passes < max_passes) {
      
      double g = 0.5*accu((x - beta) % reshape(hessianProduct(x - beta), size(x)));
      double h = dot(sort(abs(vectorise(x.tail_rows(p_rows))),
                          "descending"), lambda);
      double f = g + h;

      mat grad = reshape(hessianProduct(x - beta), size(x));

      x_tilde_old = x_tilde;

      double g_old = g;
      double t_old = t;

      // Backtracking line search
      while (true) {
        // Update coefficients
        x_tilde = x - learning_rate*grad;

        x_tilde.tail_rows(p_rows) =
          prox(x_tilde.tail_rows(p_rows), lambda*learning_rate);

        vec d = vectorise(x_tilde - x);

        g = 0.5*accu((x_tilde - beta) % reshape(hessianProduct(x_tilde - beta), size(x)));

        double q = g_old
          + dot(d, vectorise(grad))
          + (1.0/(2*learning_rate))*accu(square(d));

          if (q >= g*(1 - 1e-12)) {
            break;
          } else {
            learning_rate *= alpha;
          }
          checkUserInterrupt();
      }

      t = 0.5*(1.0 + std::sqrt(1.0 + 4.0*t_old*t_old));
      x = x_tilde + (t_old - 1.0)/t * (x_tilde - x_tilde_old);
      
      // Stop if change is x is small
      if (norm(x - x_prev, "fro") < 1e-6)
        break;
      x_prev = x;

      if (passes % 100 == 0)
        checkUserInterrupt();

      ++passes;
      
    }
    
    return x_tilde;

  }
  


};

