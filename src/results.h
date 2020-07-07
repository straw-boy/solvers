#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

struct Results {
  mat beta;
  uword passes;
  std::vector<std::vector<double>> diagnostics_loss;
  std::vector<double> time;
  double deviance;

  Results() {}

  Results(mat beta,
          uword passes,
          std::vector<std::vector<double>> diagnostics_loss,
          std::vector<double> time,
          double deviance)
    : beta(beta),
      passes(passes),
      diagnostics_loss(diagnostics_loss),
      time(time),
      deviance(deviance) {}
};
