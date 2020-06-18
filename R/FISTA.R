#' FISTA

FISTA <- function(x,
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
                  solver = c("fista", "admm"),
                  max_passes = 1e6,
                  tol_abs = 1e-5,
                  tol_rel = 1e-4,
                  tol_rel_gap = 1e-5,
                  tol_infeas = 1e-3,
                  diagnostics = FALSE,
                  verbosity = 0,
                  sigma,
                  n_sigma,
                  lambda_min_ratio
) {

  if (!missing(sigma)) {
    warning("`sigma` argument is deprecated; please use `alpha` instead")
    alpha <- sigma
  }

  if (!missing(n_sigma)) {
    warning("`n_sigma` argument is deprecated; please use `path_length` instead")
    path_length <- n_sigma
  }

  if (!missing(lambda_min_ratio)) {
    warning("`lambda_min_ratio` is deprecated; please use `alpha_min_ratio` instead")
    alpha_min_ratio <- lambda_min_ratio
  }

  ocall <- match.call()

  family <- match.arg(family)
  solver <- match.arg(solver)
  screen_alg <- match.arg(screen_alg)

  if (solver == "admm" && family != "gaussian")
    stop("ADMM solver is only supported with `family = 'gaussian'`")

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
    tol_rel_gap >= 0,
    tol_infeas >= 0,
    tol_abs >= 0,
    tol_rel >= 0,
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
                  solver = solver,
                  tol_dev_change = tol_dev_change,
                  tol_dev_ratio = tol_dev_ratio,
                  tol_rel_gap = tol_rel_gap,
                  tol_infeas = tol_infeas,
                  tol_abs = tol_abs,
                  tol_rel = tol_rel)

  fitFISTA <- if (is_sparse) sparseFISTA else denseFISTA

  if (intercept) {
    x <- cbind(1, x)
  }

  if (alpha_type %in% c("path", "user")) {
    fit <- fitFISTA(x, y, control)
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

        fit <- fitFISTA(x, y, control)

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
      fit <- fitFISTA(x, y, control)
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

  diagnostics <- if (diagnostics) setupDiagnostics(fit) else NULL

  slope_class <- switch(family,
                        gaussian = "GaussianSLOPE",
                        binomial = "BinomialSLOPE",
                        poisson = "PoissonSLOPE",
                        multinomial = "MultinomialSLOPE")

  structure(list(coefficients = coefficients,
                 nonzeros = nonzeros,
                 lambda = lambda,
                 alpha = alpha,
                 class_names = class_names,
                 passes = fit$passes,
                 execution_times = fit$execution_timings,
                 total_time = fit$total_time,
                 unique = drop(fit$n_unique),
                 deviance_ratio = drop(fit$deviance_ratio),
                 null_deviance = fit$null_deviance,
                 family = family,
                 diagnostics = diagnostics,
                 call = ocall),
            class = c(slope_class, "SLOPE"))
}
