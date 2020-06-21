#pragma once

#include <RcppArmadillo.h>
#include "family.h"
#include "../results.h"
using namespace Rcpp;
using namespace arma;

class Gaussian : public Family {
private:
  double alpha = 1.5;

public:
  template <typename... Ts>
  Gaussian(Ts... args) : Family(std::forward<Ts>(args)...) {}

  double primal(const mat& y, const mat& lin_pred)
  {
    return 0.5*pow(norm(y - lin_pred), 2);
  }

  double dual(const mat& y, const mat& lin_pred)
  {
    using namespace std;
    return 0.5*pow(norm(y, 2), 2) - 0.5*pow(norm(lin_pred, 2), 2);
  }
  
  mat pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return lin_pred - y;
  }

  mat pseudoHessian(const mat& y, const mat& lin_pred){
    return ones<vec>(y.n_rows);
  }
  
  rowvec fitNullModel(const mat& y, const uword n_classes)
  {
    return mean(y);
  }

  std::string name()
  {
    return "gaussian";
  }

};

