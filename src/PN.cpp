#include <RcppArmadillo.h>
#include <memory>
#include "results.h"
#include "families/families.h"
#include "algorithms/proximal_newton.h"
#include "standardize.h"
#include "rescale.h"
#include "regularizationPath.h"

using namespace Rcpp;
using namespace arma;

template <typename T>
List cppPN(T& x, mat& y, const List control)
{

  wall_clock outer_timer;
  outer_timer.tic();

  using std::endl;
  using std::setw;
  using std::showpoint;

  // significant digits
  Rcout.precision(4);

  auto tol_dev_ratio = as<double>(control["tol_dev_ratio"]);
  auto tol_dev_change = as<double>(control["tol_dev_change"]);
  auto max_variables = as<uword>(control["max_variables"]);

  auto diagnostics = as<bool>(control["diagnostics"]);
  auto verbosity = as<uword>(control["verbosity"]);

  // solver arguments
  auto max_passes  = as<uword>(control["max_passes"]);
  auto tol_rel_gap = as<double>(control["tol_rel_gap"]);
  auto tol_infeas  = as<double>(control["tol_infeas"]);
  auto tol_abs     = as<double>(control["tol_abs"]);
  auto tol_rel     = as<double>(control["tol_rel"]);

  auto family_choice = as<std::string>(control["family"]);
  auto intercept     = as<bool>(control["fit_intercept"]);

  auto n = x.n_rows;
  auto p = x.n_cols;
  auto m = y.n_cols;

  auto center = as<bool>(control["center"]);
  auto scale = as<std::string>(control["scale"]);

  auto y_center = as<rowvec>(control["y_center"]);
  auto y_scale = as<rowvec>(control["y_scale"]);
  rowvec x_center(p, fill::zeros);
  rowvec x_scale(p, fill::ones);

  standardize(x, x_center, x_scale, intercept, center, scale);

  auto lambda = as<vec>(control["lambda"]);
  auto alpha  = as<vec>(control["alpha"]);
  auto lambda_type = as<std::string>(control["lambda_type"]);
  auto alpha_type = as<std::string>(control["alpha_type"]);
  auto alpha_min_ratio = as<double>(control["alpha_min_ratio"]);
  auto q = as<double>(control["q"]);
  const uword path_length = alpha.n_elem;
  double alpha_max = 0;

  regularizationPath(alpha,
                     lambda,
                     alpha_max,
                     x,
                     y,
                     x_scale,
                     y_scale,
                     lambda_type,
                     alpha_type,
                     scale,
                     alpha_min_ratio,
                     q,
                     family_choice,
                     intercept);

  auto family = setupFamily(family_choice,
                            intercept,
                            diagnostics,
                            max_passes,
                            tol_rel_gap,
                            tol_infeas,
                            tol_abs,
                            tol_rel,
                            verbosity);

  cube betas(p, m, path_length, fill::zeros);
  mat beta(p, m, fill::zeros);

  uvec n_unique(path_length);

  mat linear_predictor = x*beta;

  double null_deviance = 2*family->primal(y, linear_predictor);
  vec deviance_ratios(path_length);
  vec deviances(path_length);
  double deviance_change{0};

  mat beta_prev(p, m, fill::zeros);

  uvec passes(path_length);
  std::vector<std::vector<double>> loss;
  std::vector<std::vector<double>> iteration_timings;
  std::vector<double> execution_timings;
  

  Results res;
  uword k = 0;

  wall_clock inner_timer;

  while (k < path_length) {
    inner_timer.tic();
    res = family->fitProximalNewton(x, y, lambda*alpha(k));
    passes(k) = res.passes;
    beta = res.beta;

    if (diagnostics) {
      loss.push_back(res.primals);
      iteration_timings.push_back(res.time);
    }
    
    // store coefficients and intercept
    double deviance = res.deviance;
    double deviance_ratio = 1.0 - deviance/null_deviance;
    deviances(k) = deviance;
    deviance_ratios(k) = deviance_ratio;

    deviance_change =
      k == 0 ? 0.0 : std::abs((deviances(k-1) - deviance)/deviances(k-1));

    betas.slice(k) = beta;
    beta_prev = beta;

    uword n_coefs = accu(any(beta != 0, 1));
    n_unique(k) = unique(abs(nonzeros(beta))).eval().n_elem;

    if (n_coefs > 0 && k > 0) {
      // stop path if fractional deviance change is small
      if (deviance_change < tol_dev_change || deviance_ratio > tol_dev_ratio) {
        k++;
        break;
      }
    }

    if (n_unique(k) > max_variables)
      break;

    execution_timings.push_back(inner_timer.toc());
    k++;

    checkUserInterrupt();
  }

  betas.resize(p, m, k);
  passes.resize(k);
  alpha.resize(k);
  n_unique.resize(k);
  deviances.resize(k);
  deviance_ratios.resize(k);

  rescale(betas,
          x_center,
          x_scale,
          y_center,
          y_scale,
          intercept);

  // rescale alpha depending on standardization settings
  if (scale == "l2") {
    alpha /= sqrt(n);
  } else if (scale == "sd" || scale == "none") {
    alpha /= n;
  }

  return List::create(
    Named("betas")               = wrap(betas),
    Named("passes")              = wrap(passes),
    Named("loss")                = wrap(loss),
    Named("iteration_timings")   = wrap(iteration_timings),
    Named("execution_timings")   = wrap(execution_timings),
    Named("total_time")          = wrap(outer_timer.toc()),
    Named("n_unique")            = wrap(n_unique),
    Named("deviance_ratio")      = wrap(deviance_ratios),
    Named("null_deviance")       = wrap(null_deviance),
    Named("alpha")               = wrap(alpha),
    Named("lambda")              = wrap(lambda)
  );
}

// [[Rcpp::export]]
Rcpp::List sparsePN(arma::sp_mat x,
                       arma::mat y,
                       const Rcpp::List control)
{
  return cppPN(x, y, control);
}

// [[Rcpp::export]]
Rcpp::List densePN(arma::mat x,
                      arma::mat y,
                      const Rcpp::List control)
{
  return cppPN(x, y, control);
}
