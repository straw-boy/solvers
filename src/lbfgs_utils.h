#pragma once

#include <RcppArmadillo.h>
#include "prox.h"
#include "infeasibility.h"

using namespace Rcpp;
using namespace arma;


class LBFGS {
private:
  const uword m = 7;
  double gamma = 1.0;
  double sigma = 1.0;
  
  mat S;
  mat Y;
  mat L;
  mat R;
  mat D;

  mat Q_Hv;
  mat Q_Bv;
  mat A_Hv;
  mat A_Bv;

  // Scaled Proximal parameters
  const double tol_rel_gap = 1e-5;
  const double tol_infeas = 1e-3;
  const double tol_rel_coef_change = 1e-3;
  uword max_passes = 10000;

public:

  bool updateParams(const mat& s, const mat& y) {

    if (dot(s, y) == 0)
      return true;
    
    S.insert_cols(S.n_cols, vectorise(s));
    Y.insert_cols(Y.n_cols, vectorise(y));

    if (S.n_cols > m) {
      S.shed_col(0);
      Y.shed_col(0);
    }

    gamma = dot(s, y)/dot(y, y);
    sigma = 1/gamma;
    // sigma = dot(s, y)/dot(s, s);

    // Rcout << "--------" << endl;
    mat sTy = S.t()*Y;
    // sTy.print();
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

  mat computeHv(const mat& v) {

    if (S.n_cols == 0) {
      return gamma*v;
    }

    // Rcout << " Hv main entered" << endl;

    mat Rinv = inv(R);
    mat z = vectorise(v);
    mat Q = join_rows(S, gamma*Y);

    mat ans = Q.t()*z;

    mat tmp = Rinv.t()*(D + gamma*Y.t()*Y)*Rinv;
    tmp = join_rows(tmp, -Rinv.t());
    tmp = join_cols(tmp, join_rows(-Rinv, zeros(size(R))));

    ans = tmp*ans;
    ans = Q*ans;
    ans = gamma*z + ans;

    return reshape(ans, size(v));
  }

  mat computeHv_cached(const mat& v) {

    if (S.n_cols == 0) {
      return gamma*v;
    }

    mat ans = Q_Hv.t()*v;
    ans = A_Hv*ans;
    ans = Q_Hv*ans;
    ans = gamma*v + ans;
    return ans;

  }


  mat computeBv(const mat& v) {

    if (S.n_cols == 0) {
      return sigma*v;
    }
    // Rcout << " Bv main entered" << endl;

    mat z = vectorise(v);
    mat Q = join_rows(sigma*S, Y);

    mat ans = Q.t()*z;

    mat tmp = sigma*S.t()*S;
    tmp = join_rows(tmp, L);
    tmp = join_cols(tmp, join_rows(L.t(), -D));

    // ans = solve(tmp, ans);
    ans = inv(tmp)*ans;
    ans = Q*ans;
    ans = sigma*z - ans;
    
    return ans;
    // return reshape(ans, size(v));
  }

  mat computeBv_cached(const mat& v) {

    if (S.n_cols == 0) {
      return sigma*v;
    }

    mat ans = Q_Bv.t()*v;
    ans = inv(A_Bv)*ans;
    ans = Q_Bv*ans;
    ans = sigma*v - ans;
    return ans;

  }



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

    // line search parameters
    double alpha = 0.5;

    // FISTA parameter
    double t = 1;

    uword passes = 0;
    while (passes < max_passes) {
      
      double g = 0.5*dot(x - beta, computeBv_cached(x - beta));
      double h = dot(sort(abs(vectorise(x.tail_rows(p_rows))),
                          "descending"), lambda);
      double f = g + h;

      mat grad = computeBv(x - beta);

      x_tilde_old = x_tilde;

      double g_old = g;
      double t_old = t;

      // Backtracking line search
      while (true) {
        // Update coefficients
        x_tilde = x - learning_rate*(grad);

        x_tilde.tail_rows(p_rows) =
          prox(x_tilde.tail_rows(p_rows), lambda*learning_rate);

        vec d = vectorise(x_tilde - x);

        g = 0.5*dot(x_tilde - beta, computeBv_cached(x_tilde - beta));

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

      // if(learning_rate < 1e-8)
      //   break;

      // FISTA step
      t = 0.5*(1.0 + std::sqrt(1.0 + 4.0*t_old*t_old));
      x = x_tilde + (t_old - 1.0)/t * (x_tilde - x_tilde_old);
      
      if(norm(x-x_prev) < 1e-6)
        break;
      x_prev = x;

      if (passes % 100 == 0)
        checkUserInterrupt();

      ++passes;
      
    }
    

    return x_tilde;

  }
  


};

