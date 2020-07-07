#' Code borrowed from SLOPE package for generating datasets
#' @param n Number of data points
#' @param p Number of features
#' @param q <TODO>
#' @param n_groups <TODO>
#' @param n_targets <TODO>
#' @param density Determine sparsity of the data generated. Default set to 0.
#' @param amplitude <TODO>
#' @param alpha <TODO>
#' @param response Choice of likelohood function
#' @param rho <TODO>
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

#' FISTA
#' @param fits List of outputs of any of the algorithms (FISTA, ADMM etc)
#' @export 
mergeFits <- function(fits) {
  
  f <- data.frame(solver = character(),
                  loss = double(),
                  time = double())

  for (fit in fits) {
    g  <- data.frame(solver = rep(getLabel(fit), times = length(fit$iteration_timings[[1]])),
                     loss = getLoss(fit),
                     time = fit$iteration_timings[[1]])
    f <- rbind(f,g)
  }
  return(f)
}



