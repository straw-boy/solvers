% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{randomProblem}
\alias{randomProblem}
\title{Code borrowed from SLOPE package for generating datasets}
\usage{
randomProblem(
  n = 1000,
  p = 100,
  q = 0.2,
  n_groups = NULL,
  n_targets = if (match.arg(response) == "multinomial") 3 else 1,
  density = 0,
  amplitude = if (match.arg(response) == "poisson") 1 else 3,
  alpha = 1,
  response = c("gaussian", "binomial", "poisson", "multinomial"),
  rho = 0
)
}
\arguments{
\item{n}{Number of data points}

\item{p}{Number of features}

\item{q}{Parameter controlling lambda sequence}

\item{n_groups}{Number of groups of predictors}

\item{n_targets}{Dimension of response variable}

\item{density}{Determine sparsity of the data generated. Default set to 0.}

\item{amplitude}{Scale of coefficient in output}

\item{alpha}{Standard deviation}

\item{response}{Choice of likelohood function}

\item{rho}{Correlation}
}
\description{
Code borrowed from SLOPE package for generating datasets
}
