test_that("Proximal Newton: gaussian, n>p case", {
  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- solvers::randomProblem(n, p, response="gaussian", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), opt_algo="nr")
  pn_solvers <- PN(d$x, d$y, family="gaussian", alpha=c(1.0,0.005))
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

})

test_that("Proximal Newton: gaussian, n<p case", {

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- solvers::randomProblem(n, p, response="gaussian", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="gaussian",alpha=c(1.0,0.005),opt_algo="nr")
  pn_solvers <- PN(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

})

test_that("Proximal Newton: binomial, n>p case", {
  
  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- solvers::randomProblem(n, p, response="binomial", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005),opt_algo="nr")
  pn_solvers <- PN(d$x, d$y, family="binomial",alpha=c(1.0,0.005))
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

})


test_that("Proximal Newton: binomial, n<p case", {

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- solvers::randomProblem(n, p, response="binomial", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005),opt_algo="nr")
  pn_solvers <- PN(d$x, d$y, family="binomial",alpha=c(1.0,0.005))
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

})

test_that("Proximal Newton: poisson, n>p case", {

  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- solvers::randomProblem(n, p, response="poisson", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005),opt_algo="nr")
  pn_solvers <- PN(d$x, d$y, family="poisson",alpha=c(1.0,0.005))
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

})

test_that("Proximal Newton: poisson, n<p case", {

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- solvers::randomProblem(n, p, response="poisson", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005),opt_algo="nr")
  pn_solvers <- PN(d$x, d$y, family="poisson",alpha=c(1.0,0.005))
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

})

