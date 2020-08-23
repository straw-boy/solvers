test_that("L-BFGS: gaussian, n>p case", {
  library(SLOPE)
  set.seed(326)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="gaussian")
  
  nr <- ADMM(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="gaussian",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})

test_that("L-BFGS: gaussian, n<p case", {

  library(SLOPE)
  set.seed(221)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="gaussian")
  
  nr <- ADMM(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="gaussian",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})

test_that("L-BFGS: binomial, n>p case", {
  
  library(SLOPE)
  set.seed(541)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="binomial")
  
  nr <- ADMM(d$x, d$y, family="binomial", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})


test_that("L-BFGS: binomial, n<p case", {

  library(SLOPE)
  set.seed(771)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="binomial")
  
  nr <- ADMM(d$x, d$y, family="binomial", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})

test_that("L-BFGS: poisson, n>p case", {

  library(SLOPE)
  set.seed(751)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="poisson")
  
  nr <- ADMM(d$x, d$y, family="poisson", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)
})

test_that("L-BFGS: poisson, n<p case", {

  library(SLOPE)
  set.seed(122)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="poisson")
  
  nr <- ADMM(d$x, d$y, family="poisson", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)

})

test_that("L-BFGS: multinomial, n>p case", {

  library(SLOPE)
  set.seed(46)

  n = 100
  p = 10

  d <- randomProblem(n, p, response="multinomial", density = 0.5)
  
  nr <- ADMM(d$x, d$y, family="multinomial", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="multinomial",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)
  
})

test_that("L-BFGS: multinomial, n<p case", {

  library(SLOPE)
  set.seed(67)

  n = 10
  p = 20

  d <- randomProblem(n, p, response="multinomial", density = 0.5)
  
  nr <- ADMM(d$x, d$y, family="multinomial", alpha=c(1.0,0.005), opt_algo="nr")
  lbfgs <- ADMM(d$x, d$y, family="multinomial",alpha=c(1.0,0.005), opt_algo="lbfgs")
  expect_equivalent(coef(nr), coef(lbfgs), tol = 1e-2)
  
})

