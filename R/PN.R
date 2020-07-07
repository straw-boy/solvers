#' Proximal Newton
#' @param x the design matrix, which can be either a dense
#'   matrix of the standard *matrix* class, or a sparse matrix
#'   inheriting from [Matrix::sparseMatrix]. Data frames will
#'   be converted to matrices internally.
#' @param y the response, which for `family = "gaussian"` must be numeric; for
#'   `family = "binomial"` or `family = "multinomial"`, it can be a factor.
#' @param family model family 
#' @param intercept whether to fit an intercept
#' @param center whether to center predictors or not by their mean. Defaults
#'   to `TRUE` if `x` is dense and `FALSE` otherwise.
#' @param scale type of scaling to apply to predictors.
#'   - `"l1"` scales predictors to have L1 norms of one.
#'   - `"l2"` scales predictors to have L2 norms of one.#'
#'   - `"sd"` scales predictors to have a population standard deviation one.
#'   - `"none"` applies no scaling.
#' @param alpha scale for regularization path: either a decreasing numeric
#'   vector (possibly of length 1) or a character vector; in the latter case,
#'   the choices are:
#'   - `"path"`, which computes a regularization sequence
#'     where the first value corresponds to the intercept-only (null) model and
#'     the last to the almost-saturated model, and
#'   - `"estimate"`, which estimates a *single* `alpha`
#'     using Algorithm 5 in Bogdan et al. (2015).
#' @param path_length length of regularization path; note that the path
#'   returned may still be shorter due to the early termination criteria
#'   given by `tol_dev_change`, `tol_dev_ratio`, and `max_variables`.
#' @param lambda either a character vector indicating the method used
#'   to construct the lambda path or a numeric non-decreasing
#'   vector with length equal to the number
#'   of coefficients in the model; see section **Regularization sequences**
#'   for details.
#' @param alpha_min_ratio smallest value for `lambda` as a fraction of
#'   `lambda_max`; used in the selection of `alpha` when `alpha = "path"`
#' @param q parameter controlling the shape of the lambda sequence, with
#'   usage varying depending on the type of path used and has no effect
#'   is a custom `lambda` sequence is used.
#' @param max_passes maximum number of passes (outer iterations) for solver
#' @param diagnostics whether to save diagnostics from the solver
#'   (timings and other values depending on type of solver)
#' @param screen (currently inactive) whether to use predictor screening rules 
#'   (rules that allow some predictors to be discarded prior to fitting),
#'   which improve speed greatly when the number of predictors is larger than 
#'   the number of observations.
#' @param screen_alg (currently inactive) what type of screening algorithm to use.
#'   - `"strong"` uses the set from the strong screening rule and check
#'     against the full set
#'   - `"previous"` first fits with the previous active set, then checks
#'     against the strong set, and finally against the full set if there are
#'     no violations in the strong set
#' @param verbosity level of verbosity for displaying output from the
#'   program. Not completely developed. Use 3 just for now.
#' @param tol_dev_change the regularization path is stopped if the
#'   fractional change in deviance falls below this value; note that this is
#'   automatically set to 0 if a alpha is manually entered
#' @param tol_dev_ratio the regularization path is stopped if the
#'   deviance ratio
#'   \eqn{1 - \mathrm{deviance}/\mathrm{(null-deviance)}}{1 - deviance/(null deviance)}
#'   is above this threshold
#' @param max_variables criterion for stopping the path in terms of the
#'   maximum number of unique, nonzero coefficients in absolute value in model.
#'   For the multinomial family, this value will be multiplied internally with
#'   the number of levels of the response minus one.
#' @export
PN <- function(x,
               y,
               family = c("gaussian", "binomial", "multinomial", "poisson"),
               intercept = TRUE,
               center = !inherits(x, "sparseMatrix"),
               scale = c("l2", "l1", "sd", "none"),
               alpha = c("path", "estimate"),
               lambda = c("bh", "gaussian", "oscar"),
               alpha_min_ratio = if (NROW(x) < NCOL(x)) 1e-2 else 1e-4,
               path_length = if (alpha[1] == "estimate") 1 else 20,
               q = 0.1*min(1, NROW(x)/NCOL(x)),
               screen = TRUE,
               screen_alg = c("strong", "previous"),
               tol_dev_change = 1e-5,
               tol_dev_ratio = 0.995,
               max_variables = NROW(x),
               max_passes = 100,
               diagnostics =  FALSE,
               verbosity = 0
) {

  ocall <- match.call()

  family <- match.arg(family)
  screen_alg <- match.arg(screen_alg)

  if (is.character(scale)) {
    scale <- match.arg(scale)
  } else if (is.logical(scale) && length(scale) == 1L) {
    scale <- ifelse(scale, "l2", "none")
  } else {
    stop("`scale` must be logical or a character")
  }

  n <- NROW(x)
  p <- NCOL(x)

  stopifnot(
    is.null(alpha_min_ratio) ||
      (alpha_min_ratio > 0 && alpha_min_ratio < 1),
    max_passes > 0,
    q > 0,
    q < 1,
    length(path_length) == 1,
    path_length >= 1,
    is.null(lambda) || is.character(lambda) || is.numeric(lambda),
    is.finite(max_passes),
    is.logical(diagnostics),
    is.logical(intercept),
    is.logical(center)
  )

  fit_intercept <- intercept

  # convert sparse x to dgCMatrix class from package Matrix.
  is_sparse <- inherits(x, "sparseMatrix")

  if (NROW(y) != NROW(x))
    stop("the number of samples in `x` and `y` must match")

  if (NROW(y) == 0)
    stop("`y` is empty")

  if (NROW(x) == 0)
    stop("`x` is empty")

  if (anyNA(y) || anyNA(x))
    stop("missing values are not allowed")

  if (is_sparse) {
    x <- methods::as(x, "dgCMatrix")
  } else {
    x <- as.matrix(x)
  }

  if (is_sparse && center)
    stop("centering would destroy sparsity in `x` (predictor matrix)")

  res <- preprocessResponse(family, y)
  y <- as.matrix(res$y)
  y_center <- res$y_center
  y_scale <- res$y_scale
  class_names <- res$class_names
  m <- n_targets <- res$n_targets
  response_names <- res$response_names
  variable_names <- colnames(x)
  max_variables <- max_variables*m

  if (is.null(variable_names))
    variable_names <- paste0("V", seq_len(p))
  if (is.null(response_names))
    response_names <- paste0("y", seq_len(m))

  if (is.character(alpha)) {
    alpha <- match.arg(alpha)

    if (alpha == "path") {

      alpha_type <- "auto"
      alpha <- double(path_length)

    } else if (alpha == "estimate") {

      if (family != "gaussian")
        stop("`alpha = 'estimate'` can only be used if `family = 'gaussian'`")

      alpha_type <- "estimate"
      alpha <- NULL

      if (path_length > 1)
        warning("`path_length` ignored since `alpha = 'estimate'`")
    }
  } else {
    alpha <- as.double(alpha)
    alpha_type <- "user"

    alpha <- as.double(alpha)
    path_length <- length(alpha)

    stopifnot(path_length > 0)

    if (any(alpha < 0))
      stop("`alpha` cannot contain negative values")

    if (is.unsorted(rev(alpha)))
      stop("`alpha` must be decreasing")

    if (anyDuplicated(alpha) > 0)
      stop("all values in `alpha` must be unique")

    # do not stop path early if user requests specific alpha
    tol_dev_change <- 0
    tol_dev_ratio <- 1
    max_variables <- (NCOL(x) + intercept)*m
  }

  n_lambda <- m*p

  if (is.character(lambda)) {

    lambda_type <- match.arg(lambda)

    if (lambda_type == "bhq")
      warning("'bhq' option to argument lambda has been depracted and will",
              "will be defunct in the next release; please use 'bh' instead")

    lambda <- double(n_lambda)

  } else {

    lambda_type <- "user"
    lambda <- as.double(lambda)

    if (length(lambda) != m*p)
      stop("`lambda` must be as long as there are variables")

    if (is.unsorted(rev(lambda)))
      stop("`lambda` must be non-increasing")

    if (any(lambda < 0))
      stop("`lambda` cannot contain negative values")
  }

  control <- list(family = family,
                  fit_intercept = fit_intercept,
                  is_sparse = is_sparse,
                  scale = scale,
                  center = center,
                  path_length = path_length,
                  n_targets = n_targets,
                  screen = screen,
                  screen_alg = screen_alg,
                  alpha = alpha,
                  alpha_type = alpha_type,
                  lambda = lambda,
                  lambda_type = lambda_type,
                  alpha_min_ratio = alpha_min_ratio,
                  q = q,
                  y_center = y_center,
                  y_scale = y_scale,
                  max_passes = max_passes,
                  diagnostics = diagnostics,
                  verbosity = verbosity,
                  max_variables = max_variables,
                  tol_dev_change = tol_dev_change,
                  tol_dev_ratio = tol_dev_ratio,
                  tol_rel_gap = 1e-5,
                  tol_infeas = 1e-3,
                  tol_abs = 1e-5,
                  tol_rel = 1e-4)

  fitPN <- if (is_sparse) sparsePN else densePN

  if (intercept) {
    x <- cbind(1, x)
  }

  if (alpha_type %in% c("path", "user")) {
    fit <- fitPN(x, y, control)
  } else {
    # estimate the noise level, if possible
    if (is.null(alpha) && n >= p + 30)
      alpha <- estimateNoise(x, y)

    # run the solver, iteratively if necessary.
    if (is.null(alpha)) {
      # Run Algorithm 5 of Section 3.2.3. (Bogdan et al.)
      if (intercept)
        selected <- 1
      else
        selected <- integer(0)

      repeat {
        selected_prev <- selected

        alpha <- estimateNoise(x[, selected, drop = FALSE], y, intercept)
        control$alpha <- alpha

        fit <- fitPN(x, y, control)

        selected <- which(abs(drop(fit$betas)) > 0)

        if (fit_intercept)
          selected <- union(1, selected)

        if (identical(selected, selected_prev))
          break

        if (length(selected) + 1 >= n)
          stop("selected >= n-1 variables; cannot estimate variance")
      }
    } else {
      control$alpha <- alpha
      fit <- fitPN(x, y, control)
    }
  }

  lambda <- fit$lambda
  alpha <- fit$alpha
  path_length <- length(alpha)
  beta <- fit$betas
  nonzeros <- apply(beta, c(2, 3), function(x) abs(x) > 0)
  coefficients <- beta

  if (fit_intercept) {
    nonzeros <- nonzeros[-1, , , drop = FALSE]
    dimnames(coefficients) <- list(c("(Intercept)", variable_names),
                                   response_names[1:n_targets],
                                   paste0("p", seq_len(path_length)))
  } else {
    dimnames(coefficients) <- list(variable_names,
                                   response_names[1:n_targets],
                                   paste0("p", seq_len(path_length)))
  }

  slope_class <- switch(family,
                        gaussian = "GaussianSLOPE",
                        binomial = "BinomialSLOPE",
                        poisson = "PoissonSLOPE",
                        multinomial = "MultinomialSLOPE")

  structure(list(solver = "PN",
                 coefficients = coefficients,
                 nonzeros = nonzeros,
                 lambda = lambda,
                 alpha = alpha,
                 class_names = class_names,
                 loss = fit$loss,
                 iteration_timings = fit$iteration_timings,
                 passes = fit$passes,
                 execution_times = fit$execution_times,
                 total_time = fit$total_time,
                 unique = drop(fit$n_unique),
                 deviance_ratio = drop(fit$deviance_ratio),
                 null_deviance = fit$null_deviance,
                 family = family,
                 call = ocall),
            class = c(slope_class, "SLOPE"))
}
