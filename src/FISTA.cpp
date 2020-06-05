#include <RcppArmadillo.h>
#include <memory>
#include "results.h"
#include "families/families.h"
#include "standardize.h"
#include "rescale.h"
#include "regularizationPath.h"

using namespace Rcpp;
using namespace arma;

template <typename T>
List cppFISTA(T& x, mat& y, const List control)
{
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
  auto solver      = as<std::string>(control["solver"]);
  auto max_passes  = as<uword>(control["max_passes"]);
  auto tol_rel_gap = as<double>(control["tol_rel_gap"]);
  auto tol_infeas  = as<double>(control["tol_infeas"]);
  auto tol_abs     = as<double>(control["tol_abs"]);
  auto tol_rel     = as<double>(control["tol_rel"]);

  auto family_choice = as<std::string>(control["family"]);
  auto intercept     = as<bool>(control["fit_intercept"]);
  auto screen        = as<bool>(control["screen"]);
  auto screen_alg    = as<std::string>(control["screen_alg"]);

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

  uword n_variables = 0;
  uvec n_unique(path_length);

  mat linear_predictor = x*beta;

  double null_deviance = 2*family->primal(y, linear_predictor);
  vec deviance_ratios(path_length);
  vec deviances(path_length);
  double deviance_change{0};

  mat beta_prev(p, m, fill::zeros);

  uvec passes(path_length);
  std::vector<std::vector<double>> primals;
  std::vector<std::vector<double>> duals;
  std::vector<std::vector<double>> timings;
  std::vector<unsigned> violations;
  std::vector<std::vector<unsigned>> violation_list;

  mat linear_predictor_prev(n, m);
  mat gradient_prev(p, m);
  mat pseudo_gradient_prev(n, m);
  mat L(n,m);
  mat U(n,m);
  bool factorized = false;

  Results res;

  uword k = 0;

  vec z(p,fill::zeros);
  vec u(p,fill::zeros);
  vec xTy(p,fill::zeros);
  double rho=1.0;
  while (k < path_length) {

    res = family->fit(x, y, beta, z, u, L, U, xTy, lambda*alpha(k), rho, solver);
    passes(k) = res.passes;
    beta = res.beta;

    if (diagnostics) {
      primals.push_back(res.primals);
      duals.push_back(res.duals);
      timings.push_back(res.time);
      violation_list.push_back(violations);
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
    n_variables = n_coefs;
    n_unique(k) = unique(abs(nonzeros(beta))).eval().n_elem;

    if (verbosity >= 1)
      Rcout << showpoint
            << "penalty: "      << setw(2) << k
            << ", dev: "        << setw(7) << deviance
            << ", dev ratio: "  << setw(7) << deviance_ratio
            << ", dev change: " << setw(7) << deviance_change
            << ", n var: "      << setw(5) << n_variables
            << ", n unique: "   << setw(5) << n_unique(k)
            << endl;

    if (n_coefs > 0 && k > 0) {
      // stop path if fractional deviance change is small
      if (deviance_change < tol_dev_change || deviance_ratio > tol_dev_ratio) {
        k++;
        break;
      }
    }

    if (n_unique(k) > max_variables)
      break;

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
    Named("primals")             = wrap(primals),
    Named("duals")               = wrap(duals),
    Named("time")                = wrap(timings),
    Named("n_unique")            = wrap(n_unique),
    Named("deviance_ratio")      = wrap(deviance_ratios),
    Named("null_deviance")       = wrap(null_deviance),
    Named("alpha")               = wrap(alpha),
    Named("lambda")              = wrap(lambda)
  );
}

// [[Rcpp::export]]
Rcpp::List sparseFISTA(arma::sp_mat x,
                       arma::mat y,
                       const Rcpp::List control)
{
  return cppFISTA(x, y, control);
}

// [[Rcpp::export]]
Rcpp::List denseFISTA(arma::mat x,
                      arma::mat y,
                      const Rcpp::List control)
{
  return cppFISTA(x, y, control);
}
