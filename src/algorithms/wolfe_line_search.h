#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;

// Step length lies in [t_low, t_high]. This function reduces the interval 
// until a satisfactory estimate is found. This algorithm maintains some invariants.
// Refer to the following book for more details.
// Reference - Pg.60 of Numerical Optimization, Nocedal & Wright.
template <typename T>
double Family::zoom(const T& x,
                    const mat& y,
                    const double rho,
                    const mat& u, 
                    const mat& z,
                    const mat& d,
                    double t_low,
                    double t_high)
{
  
  // Objective, gradient and decrement for step length = 0
  mat lin_pred = x*z;
  const double f0 = primal(y, lin_pred) + 0.5*rho*pow(norm(z-u), 2);
  const mat g0 = gradient(x, y, lin_pred) + rho*(z-u);
  const double dec0 = dot(g0, d);

  const double c1 = 1e-4;
  const double c2 = 0.9;

  uword max_iter = 20;
  uword iter = 0;

  double t_mid;

  while (iter < max_iter) {

    // Objective for step length = t_low
    double f_low = primal(y, x*(z+t_low*d)) + 
                   0.5*rho*pow(norm(z+t_low*d-u), 2);

    t_mid = 0.5*(t_low+t_high);

    mat z_mid = z+t_mid*d;

    // Objective for step length = t_mid
    mat lin_pred = x*z_mid;
    double f_mid = primal(y, lin_pred) + 
                   0.5*rho*pow(norm(z_mid-u), 2);

    // t_mid violates the sufficient decrease condition, reduce upper bound.
    if ((f_mid > f0 + c1*t_mid*dec0) || (f_mid >= f_low)) {
      t_high = t_mid;
    } else {

      // Decrement for step length = t_mid
      double dec_mid = dot(gradient(x, y, lin_pred) + rho*(z_mid-u), d);

      // t_mid satisfies the curvature condition
      if (std::abs(dec_mid) <= -c2*dec0) {
        return t_mid;
      }

      // t_high violates one of the invariants
      if (dec_mid*(t_high-t_low) >= 0) {
        t_high = t_low;
      }

      // t_mid has lower value objective value than t_low
      // since control entered 'else' block
      t_low = t_mid;
    }

    iter++;
  }

  return t_mid;

}

// Line search algorithm to find step length satisfying Wolfe conditions. 
// Reference - Pg.59 of Numerical Optimization, Nocedal & Wright.
template <typename T>
double Family::wolfeLineSearch(const T& x,
                               const mat&y,
                               const double rho,
                               const mat&u,
                               const mat& z,
                               const mat& d)
{

  mat lin_pred = x*z;

  // Objective, gradient and decrement for step length = 0
  const double f0 = primal(y, lin_pred) + 0.5*rho*pow(norm(z-u), 2);
  const mat g0 = gradient(x, y, lin_pred) + rho*(z-u);
  const double dec0 = dot(g0, d);


  double t_prev = 0;
  double t = 1;
  const double t_max = 2;

  // Wolfe conditions parameters
  const double c1 = 1e-4;
  const double c2 = 0.9;

  double f = f0;

  uword max_iter = 30;
  uword iter = 1;

  while (iter < max_iter) {

    mat z_new = z + t*d;

    lin_pred = x*z_new;
    double f_new = primal(y, lin_pred) + 0.5*rho*pow(norm(z_new-u), 2);

    // Current t violates sufficient decrease condtion,
    // answer lies in interval [t_prev,t] 
    if ((f_new > f0 + c1*t*dec0) || ((f_new >= f) && (iter > 1))) {
      return zoom(x, y, rho, u, z, d, t_prev, t);
    }

    mat g_new = gradient(x, y, lin_pred) + rho*(z_new-u);
    double dec_new = dot(g_new, d);

    // t satisfies the curvature condition
    if (std::abs(dec_new) <= -c2*dec0) {
      return t;
    }

    // t violates curvature condition (but satisfies sufficient decrease condition)
    // answer lies in interval [t_prev,t] 
    if (dec_new >= 0) {
      return zoom(x, y, rho, u, z, d, t, t_prev);
    }

    t_prev = t;
    f = f_new;

    // Setting t to any value in (t,t_max)
    t = t + (t_max-t)*0.8;

    iter++;

  }

  return t;

}
