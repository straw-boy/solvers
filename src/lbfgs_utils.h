#pragma once

#include <RcppArmadillo.h>
#include "prox.h"

using namespace Rcpp;
using namespace arma;


class LBFGS {
private:
  const uword k = 10;
  const vec lambda;
  double gamma = 1.0;
  double sigma = 1.0;
  mat S;
  mat Y;
  mat L;
  mat R;
  mat D;

  // Scaled Proximal parameters
  const double alpha = 1.5;
  const double rho = 1.0;
  const double tol_abs = 1e-6;
  const double tol_rel = 1e-5;
  uword max_iter = 500000;

public:
  LBFGS(const vec lambda) : lambda(lambda) {}

  bool updateParams(const mat& s, const mat& y) {

    if (dot(s, y) == 0)
      return true;
    
    S.insert_cols(S.n_cols, vectorise(s));
    Y.insert_cols(Y.n_cols, vectorise(y));

    if (S.n_cols > k) {
      S.shed_col(0);
      Y.shed_col(0);
    }

    gamma = dot(s, y)/dot(y, y);
    sigma = dot(s, y)/dot(s, s);

    mat sTy = S.t()*Y;
    L = trimatl(sTy, -1);
    R = trimatu(sTy);
    D = diagmat(sTy);

    return false;

  }

  mat computeHv(const mat& v) {

    if (S.n_cols == 0) {
      return gamma*v;
    }

    Rcout << " Hv main entered" << endl;

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

  mat computeBv(const mat& v) {

    if (S.n_cols == 0) {
      return sigma*v;
    }
    Rcout << " Bv main entered" << endl;

    mat z = vectorise(v);
    mat Q = join_rows(Y, sigma*S);

    mat ans = Q.t()*z;

    mat tmp = sigma*S.t()*S;
    tmp = join_rows(tmp, L);
    tmp = join_cols(tmp, join_rows(L.t(), -D));

    ans = solve(tmp, ans);
    ans = Q*ans;
    ans = sigma*z - ans;
    
    return reshape(ans, size(v));
  }

  mat scaled_prox(const mat& beta) {
    uword p = beta.n_rows;
    uword m = beta.n_cols;
    uword pmi = lambda.n_elem;
    uword p_rows = pmi/m;

    mat x(beta);
    mat z(x);
    mat u(x);

    uword iter = 0;

    mat Bbeta = computeBv(beta);

    lambda.print();

    while (iter < max_iter) {
      iter++;
      Rcout << "pass: " << iter << " SP obj: "
            << 0.5*dot(x-beta,computeBv(x-beta)) +
              dot(sort(abs(vectorise(beta.tail_rows(p_rows))),
                          "descending"), lambda) << endl; 

      x = beta + z - u + Bbeta/rho + rho*computeHv(z - u);
      
      mat z_old = z;
      mat x_hat = alpha*x + (1 - alpha)*z_old;

      z = x_hat + u;

      z.tail_rows(p_rows) = prox(z.tail_rows(p_rows), lambda/rho);

      u += (x_hat-z);

      double r_norm = norm(x - z);
      double s_norm = norm(rho*(z - z_old));

      double eps_primal = std::sqrt(p)*tol_abs + tol_rel*std::max(norm(x), norm(z));
      double eps_dual = std::sqrt(p)*tol_abs + tol_rel*norm(rho*u);

      if (r_norm < eps_primal && s_norm < eps_dual)
          break;

      if (iter % 1000 == 0)
        Rcpp::checkUserInterrupt();
    }
    
    return x;
  }


};