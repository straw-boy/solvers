#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

struct Results {
  mat beta;
  uword passes;
  std::vector<std::vector<double>> diagnosticsLoss;
  std::vector<double> time;
  double deviance;

  Results() {}

  Results(mat beta,
          uword passes,
          std::vector<std::vector<double>> diagnosticsLoss,
          std::vector<double> time,
          double deviance)
    : beta(beta),
      passes(passes),
      diagnosticsLoss(diagnosticsLoss),
      time(time),
      deviance(deviance) {}
};
