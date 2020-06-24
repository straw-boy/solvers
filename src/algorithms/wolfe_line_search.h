#pragma once

#include <RcppArmadillo.h>
#include "../families/family.h"

using namespace Rcpp;
using namespace arma;



template <typename T>
double Family::zoom(const T& x, const mat&y, const double rho, const mat&u, 
                    const mat& z, const mat& d, double t_low, double t_high){
  
  mat lin_pred = x*z;
  const double f0 = primal(y,lin_pred)+0.5*rho*pow(norm(z-u),2);
  const mat g0 = gradient(x,y,lin_pred)+rho*(z-u);
  const double dec0 = dot(g0,d);

  const double c1 = 1e-4;
  const double c2 = 0.9;

  uword max_iter=30;
  uword iter=0;

  while(true){

    // Rcout << t_low << " " << t_high << endl;

    double f_low = primal(y,x*(z+t_low*d)) + 
                   0.5*rho*pow(norm(z+t_low*d-u),2);

    double t_mid = 0.5*(t_low+t_high);
    mat z_mid = z+t_mid*d;
    mat lin_pred = x*z_mid;
    double f_mid = primal(y,lin_pred) + 
                   0.5*rho*pow(norm(z_mid-u),2);

    if((f_mid > f0 + c1*t_mid*dec0) || (f_mid >= f_low)){
      t_high = t_mid;
    }
    else{
      double dec_mid = dot(gradient(x,y,lin_pred)+rho*(z_mid-u),d);

      if (std::abs(dec_mid) <= -c2*dec0){
        return t_mid;
      }
      if(dec_mid*(t_high-t_low)>=0){
        t_high = t_low;
      }
      t_low = t_mid;
    }
    iter++;
    if(iter>max_iter){
      // Rcout << "Max iter of zoom ";
      return t_mid;
    }
  }
  

}





template <typename T>
double Family::wolfe_line_search(const T& x, const mat&y, const double rho,
                                 const mat&u, const mat& z, const mat& d)
{

  mat lin_pred = x*z;
  const double f0 = primal(y,lin_pred)+0.5*rho*pow(norm(z-u),2);
  const mat g0 = gradient(x,y,lin_pred)+rho*(z-u);
  const double dec0 = dot(g0,d);


  double t_prev = 0;
  double t = 1;
  const double t_max = 20;
  const double c1 = 1e-4;
  const double c2 = 0.9;

  const double f = f0;

  uword max_iter=30;
  uword iter=1;

  while(true){
    mat z_new = z + t*d;

    lin_pred = x*z_new;
    double f_new = primal(y,lin_pred)+0.5*rho*pow(norm(z_new-u),2);

    if((f_new > f0 + c1*t*dec0) || ((f_new >= f) && (iter>1))){
      // Rcout << "wolfe case 1" << endl;
      return zoom(x,y,rho,u,z,d,t_prev,t);
    }

    mat g_new = gradient(x,y,lin_pred)+rho*(z_new-u);
    double dec_new = dot(g_new,d);

    if(std::abs(dec_new) <= -c2*dec0){
      // Rcout << "wolfe case 2" << endl;
      return t;
    }

    if(dec_new>=0){
      // Rcout << "wolfe case 3" << endl;
      return zoom(x,y,rho,u,z,d,t,t_prev);
    }
    t_prev = t;
    
    if(iter > max_iter){
      // Rcout << "Max iter of wolfe" << endl;
      return t;
    }
    t = t + (t_max-t)*0.8;
    Rcout << "t is now " << t << endl;

    iter++;

  }


}
