#' Code borrowed from SLOPE package for generating datasets
#' @param n Number of data points
#' @param p Number of features
#' @param q Parameter controlling lambda sequence
#' @param n_groups Number of groups of predictors
#' @param n_targets Dimension of response variable
#' @param density Determine sparsity of the data generated. Default set to 0.
#' @param amplitude Scale of coefficient in output
#' @param alpha Standard deviation
#' @param response Choice of likelohood function
#' @param rho Correlation
#' @export
randomProblem <-
  function(n = 1000,
           p = 100,
           q = 0.2,
           n_groups = NULL,
           n_targets = if (match.arg(response) == "multinomial") 3 else 1,
           density = 0,
           amplitude = if (match.arg(response) == "poisson") 1 else 3,
           alpha = 1,
           response = c("gaussian", "binomial", "poisson", "multinomial"),
           rho = 0) {
  m <- n_targets

  if (density == 1) {
    x <- matrix(stats::rnorm(n*p), n)
  } else {
    x <- Matrix::rsparsematrix(n, p, density)
  }

  if (rho > 0)
    x <- sqrt(1-rho) * x + sqrt(rho) * matrix(stats::rnorm(n), n, p)

  if (!is.null(n_groups)) {
    groups <- rep(seq_len(n_groups), each = ceiling(m*p/n_groups),
                  length.out = p*m)
    nonzero <- which(groups %in% seq_len(max(floor(n_groups*q), 1)))
  } else {
    groups <- NA
    nonzero <- sample(p*m, max(floor(q*p*m), 1))
  }

  signs <- sample(c(-1, 1), p*m, replace = TRUE)

  beta <- signs * amplitude * (1:(p*m) %in% nonzero)

  y <- switch(match.arg(response),
              gaussian = x %*% beta + stats::rnorm(n, sd = alpha),
              binomial = {
                y <- x %*% beta + stats::rnorm(n, sd = alpha)
                (sign(y) + 1)/2
              },
              multinomial = {
                beta <- matrix(beta, p, m)
                lin_pred_exp <- exp(x %*% beta)
                prob <- lin_pred_exp/rowSums(lin_pred_exp)
                y <- apply(prob, 1, function(x) sample(1:m, 1, prob = x))
              },
              poisson = {
                lambda <- as.double(exp(x %*% beta))
                y <- stats::rpois(n, lambda)
              })

  dimnames(x) <- list(seq_len(nrow(x)),
                      paste0("V", seq_len(ncol(x))))

  list(x = x,
       y = as.double(y),
       beta = beta,
       groups = groups,
       nonzero = nonzero,
       q = q)
}

getLabel <- function(fit) {
  if (fit$solver == "ADMM")
    paste("ADMM(",toupper(fit$opt_algo),")",sep='')
  else fit$solver
} 

getLoss <- function(fit) {
  if (fit$solver == "FISTA")
    fit$primal[[1]]
  else fit$loss[[1]]
} 

#' Merges all the fits and returns a single dataframe containing 
#' the base 10 log of loss, time and iteration of all the solver fits provided. 
#' @param fits List of outputs of any of the algorithms (FISTA, ADMM etc)
#' @param cutoff_time If any of the solvers run for time longer than 'cutoff_time',
#'        their output is considered only until 'cutoff_time'. This is done to avoid
#'        skewed plots due to one of the solvers taking huge time. Default value is 10000s
#' @export 
mergeFits <- function(fits, cutoff_time = 10000) 
{ 
  f <- data.frame(solver = character(),
                  loss = double(),
                  time = double(),
                  iterations = integer())

  for (fit in fits) {
    cutoff_idx <- findInterval(cutoff_time, fit$iteration_timings[[1]])
    cutoff_idx <- max(cutoff_idx, 1)
    g  <- data.frame(solver = rep(getLabel(fit), times = cutoff_idx),
                    loss = log10(getLoss(fit)[1:cutoff_idx]),
                    time = ((fit$iteration_timings[[1]])[1:cutoff_idx]),
                    iterations = seq(1, cutoff_idx))
    f <- rbind(f,g)
  }
  return(f)
}

#' Runs all the solvers on (x, y) training data with SLOPE parameter alpha,
#' prints the total time in each case, and returns the merged data frame.
#' It is important to note that ADMM(BFGS) is not run when choice of family is poisson
#' due to exploding gradients in large scale problems.
#' @param x the design matrix, which can be either a dense
#'   matrix of the standard *matrix* class, or a sparse matrix
#'   inheriting from [Matrix::sparseMatrix]. Data frames will
#'   be converted to matrices internally.
#' @param y the response, which for `family = "gaussian"` must be numeric; for
#'   `family = "binomial"` or `family = "multinomial"`, it can be a factor.
#' @param family model family 
#' @param alpha parameter used for SLOPE regularization
#' @param path_length The regularization path length. By default, it is one.
#'        'alpha' is ignored if path_length is not 1, in which case, list of alpha 
#'        and the length of the list is returned. 
#' @export 
getBenchmarks <- function(x,
                          y,
                          family = c("gaussian", "binomial", "multinomial", "poisson"),
                          alpha = 0.01,
                          path_length = 1) 
{ 
  # If path_length is not 1, no need to store diagnostic values 
  # as they won't be meaningful for plotting
  if (path_length != 1) {
    min_path_length <- 100
    fista_fit <- FISTA(x, y, family=family, path_length=path_length)
    print(paste("Time taken by FISTA        : ", fista_fit$total_time))
    min_path_length <- min(min_path_length, length(fista_fit$alpha))

    admm_nr_fit <- ADMM(x, y, family=family, opt_algo="nr", path_length=path_length)
    print(paste("Time taken by ADMM(NR)     : ", admm_nr_fit$total_time))
    min_path_length <- min(min_path_length, length(admm_nr_fit$alpha))
    
    admm_bfgs_fit <- ADMM(x, y, family=family, opt_algo="bfgs",  path_length=path_length)
    print(paste("Time taken by ADMM(BFGS)   : ", admm_bfgs_fit$total_time))
    min_path_length <- min(min_path_length, length(admm_bfgs_fit$alpha))
    
    admm_lbfgs_fit <- ADMM(x, y, family=family, opt_algo="lbfgs", path_length=path_length)
    print(paste("Time taken by ADMM(L-BFGS) : ", admm_lbfgs_fit$total_time))
    min_path_length <- min(min_path_length, length(admm_lbfgs_fit$alpha))
    
    pn_fit <- PN(x, y, family=family, hessian_calc="exact", path_length=path_length)
    print(paste("Time taken by PN           : ", pn_fit$total_time))
    min_path_length <- min(min_path_length, length(pn_fit$alpha))

    return(list(alpha = fista_fit$alpha,
                path_length = min_path_length))
  }
  fista_fit <- FISTA(x, y, family=family, alpha=alpha, diagnostics=TRUE)
  admm_nr_fit <- ADMM(x, y, family=family, opt_algo="nr", alpha=alpha, diagnostics=TRUE)
  admm_bfgs_fit <- ADMM(x, y, family=family, opt_algo="bfgs", alpha=alpha, diagnostics=TRUE)
  admm_lbfgs_fit <- ADMM(x, y, family=family, opt_algo="lbfgs", alpha=alpha, diagnostics=TRUE)
  pn_fit <- PN(x, y, family=family, alpha=alpha, hessian_calc="exact", diagnostics=TRUE)

  fits <- list(fista_fit, admm_nr_fit, admm_bfgs_fit, admm_lbfgs_fit, pn_fit)

  # Finding out the median total_time
  solver_timings <- c()
  for (fit in fits) {
    solver_timings = cbind(solver_timings, fit$total_time)
  }
  sort(solver_timings)

  print(paste("Time taken by FISTA        : ", fista_fit$total_time))
  print(paste("Time taken by ADMM(NR)     : ", admm_nr_fit$total_time))
  print(paste("Time taken by ADMM(BFGS)   : ", admm_bfgs_fit$total_time))
  print(paste("Time taken by ADMM(L-BFGS) : ", admm_lbfgs_fit$total_time))
  print(paste("Time taken by PN           : ", pn_fit$total_time))

  f <- mergeFits(fits, cutoff_time = solver_timings[2])
  return(f)
}

