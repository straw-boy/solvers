#pragma once

#include <RcppArmadillo.h>
#include "prox.h"
#include "infeasibility.h"

using namespace Rcpp;
using namespace arma;


class LBFGS {
private:
  const uword m = 15;
  const vec lambda;
  double gamma = 1.0;
  double sigma = 1.0;
  mat S;
  mat Y;
  mat L;
  mat R;
  mat D;

  // Scaled Proximal parameters
  const double tol_rel_gap = 1e-5;
  const double tol_infeas = 1e-3;
  const double tol_rel_coef_change = 1e-3;
  uword max_passes = 5000;

public:
  LBFGS(const vec lambda) : lambda(lambda) {}

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

  // mat computeH(int p) {
  //   mat I(p, p, fill::zeros);
  //   I.diag() += 1;

  //   if (S.n_cols == 0) {
  //     return gamma*I;
  //   }

  //   mat Q = join_rows(S, gamma*Y);
  //   mat ans = Q.t();

  //   mat Rinv = inv(R);
  //   mat tmp = Rinv.t()*(D + gamma*Y.t()*Y)*Rinv;
  //   tmp = join_rows(tmp, -Rinv.t());
  //   tmp = join_cols(tmp, join_rows(-Rinv, zeros(size(R))));

  //   tmp = tmp*ans;
  //   tmp = Q*tmp;

  //   return gamma*I + tmp;
    
  // }

  // mat computeHv2(const mat& v) {

  //   int l = m;

  //   mat q = v;

  //   double alpha[l];
  //   for (int j = std::min((int)S.n_cols, l)-1; j >= 0; j--) {
  //     alpha[j] = (1/dot(S.col(j),Y.col(j)))*dot(S.col(j), q);
  //     q = q - alpha[j]*Y.col(j);
  //   }

  //   q = gamma*q;

  //   for (int j = 0; j < std::min((int)S.n_cols, l); j++) {
  //     q += S.col(j)*(alpha[j]-(1/dot(S.col(j),Y.col(j)))*dot(Y.col(j), q));
  //   }
  //   // Two Loop recursion ends

  //   return q;

    
  // }


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



  mat scaled_prox(const mat& beta) {
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
      
      double g = 0.5*dot(x - beta, computeBv(x - beta));
      double h = dot(sort(abs(vectorise(x.tail_rows(p_rows))),
                          "descending"), lambda);
      double f = g + h;
      double G = 0.5*pow(norm(beta, 2), 2) - 0.5*pow(norm(x, 2), 2);

      // if (passes % 100 == 0)
      //   Rcout << " SP: pass: " << passes << " obj: " << g << " + " << h  << " = " << g+h << endl;

      mat grad = computeBv(x - beta);

      double infeas =
        lambda.n_elem > 0.0 ? infeasibility(grad.tail_rows(p_rows), lambda) : 0.0;

      double small = std::sqrt(datum::eps);

      bool optimal =
        (std::abs(f - G)/std::max(small, std::abs(f)) < tol_rel_gap);

      bool feasible =
        lambda.n_elem > 0.0 ? infeas <= std::max(small, tol_infeas*lambda(0))
                            : true;

      // check change in coefficients
      double max_change = abs(vectorise(x - x_prev)).max();
      double max_size = abs(vectorise(beta)).max();

      bool all_zero  =
        (max_size == 0.0) && (max_change == 0.0);
      bool small_change =
        (max_size != 0.0) && (max_change/max_size <= tol_rel_coef_change);

      bool small_coef_change = all_zero || small_change;
      
      if (optimal && feasible && passes > 0 && small_coef_change)
        break;
        
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

        g = 0.5*dot(x_tilde - beta, computeBv(x_tilde - beta));

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
      x_prev = x;

      if (passes % 100 == 0)
        checkUserInterrupt();

      ++passes;
      
    }

    if (p_rows != p) {
      Rcout << " sp beta diff : " << x_tilde(0) - beta(0) << endl;
      // x_tilde(0) = beta(0);
      // x(0) = beta(0);
    }
    Rcout << "SP passes : " << passes  << endl;

    return x_tilde;

  }
  


};