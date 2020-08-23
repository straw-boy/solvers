test_that("BFGS: gaussian, n>p case", {
  library(SLOPE)
  set.seed(754)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="gaussian")
  
  nr <- ADMM(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), opt_algo="nr")
  bfgs <- ADMM(d$x, d$y, family="gaussian",alpha=c(1.0,0.005), opt_algo="bfgs")
  expect_equivalent(coef(nr), coef(bfgs), tol = 1e-2)

})

test_that("BFGS: gaussian, n<p case", {

  library(SLOPE)
  set.seed(72)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="gaussian")
  
  nr <- ADMM(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), opt_algo="nr")
  bfgs <- ADMM(d$x, d$y, family="gaussian",alpha=c(1.0,0.005), opt_algo="bfgs")
  expect_equivalent(coef(nr), coef(bfgs), tol = 1e-2)

})

test_that("BFGS: binomial, n>p case", {
  
  library(SLOPE)
  set.seed(25)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="binomial")
  
  nr <- ADMM(d$x, d$y, family="binomial", alpha=c(1.0,0.005), opt_algo="nr")
  bfgs <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005), opt_algo="bfgs")
  expect_equivalent(coef(nr), coef(bfgs), tol = 1e-2)

})


test_that("BFGS: binomial, n<p case", {

  library(SLOPE)
  set.seed(157)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="binomial")
  
  nr <- ADMM(d$x, d$y, family="binomial", alpha=c(1.0,0.005), opt_algo="nr")
  bfgs <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005), opt_algo="bfgs")
  expect_equivalent(coef(nr), coef(bfgs), tol = 1e-2)

})

test_that("BFGS: poisson, n>p case", {

  library(SLOPE)
  set.seed(122)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="poisson")
  
  nr <- ADMM(d$x, d$y, family="poisson", alpha=c(1.0,0.005), opt_algo="nr")
  bfgs <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005), opt_algo="bfgs")
  expect_equivalent(coef(nr), coef(bfgs), tol = 1e-2)
})

test_that("BFGS: poisson, n<p case", {

  library(SLOPE)
  set.seed(351)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="poisson")
  
  nr <- ADMM(d$x, d$y, family="poisson", alpha=c(1.0,0.005), opt_algo="nr")
  bfgs <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005), opt_algo="bfgs")
  expect_equivalent(coef(nr), coef(bfgs), tol = 1e-2)

})

test_that("BFGS: multinomial, n>p case", {

  library(SLOPE)
  set.seed(5)

  n = 100
  p = 10

  d <- randomProblem(n, p, response="multinomial", density = 0.5)
  
  nr <- ADMM(d$x, d$y, family="multinomial", alpha=c(1.0,0.005), opt_algo="nr")
  bfgs <- ADMM(d$x, d$y, family="multinomial",alpha=c(1.0,0.005), opt_algo="bfgs")
  expect_equivalent(coef(nr), coef(bfgs), tol = 1e-2)
  
})

test_that("BFGS: multinomial, n<p case", {

  library(SLOPE)
  set.seed(331)

  n = 10
  p = 20

  d <- randomProblem(n, p, response="multinomial", density = 0.5)
  
  nr <- ADMM(d$x, d$y, family="multinomial", alpha=c(1.0,0.005), opt_algo="nr")
  bfgs <- ADMM(d$x, d$y, family="multinomial",alpha=c(1.0,0.005), opt_algo="bfgs")
  expect_equivalent(coef(nr), coef(bfgs), tol = 1e-2)
  
})

