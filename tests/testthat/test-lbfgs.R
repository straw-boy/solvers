test_that("L-BFGS: gaussian, n>p case", {
  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="gaussian")
  
  nr <- ADMM(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="gaussian",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})

test_that("L-BFGS: gaussian, n<p case", {

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="gaussian")
  
  nr <- ADMM(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="gaussian",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})

test_that("L-BFGS: binomial, n>p case", {
  
  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="binomial")
  
  nr <- ADMM(d$x, d$y, family="binomial", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})


test_that("L-BFGS: binomial, n<p case", {

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="binomial")
  
  nr <- ADMM(d$x, d$y, family="binomial", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})

test_that("L-BFGS: poisson, n>p case", {

  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="poisson")
  
  nr <- ADMM(d$x, d$y, family="poisson", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)
})

test_that("L-BFGS: poisson, n<p case", {

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="poisson")
  
  nr <- ADMM(d$x, d$y, family="poisson", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})

