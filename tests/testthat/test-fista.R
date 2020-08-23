test_that("FISTA: gaussian, n>p case", {
  library(SLOPE)
  set.seed(131)

  n = 100
  p = 10

  d <- randomProblem(n, p, response="gaussian", density = 0.5)
  
  fista_solvers <- FISTA(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: gaussian, n<p case", {

  library(SLOPE)
  set.seed(541)

  n = 5
  p = 20

  d <- randomProblem(n, p, response="gaussian", density = 0.5)
  
  fista_solvers <- FISTA(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})


test_that("FISTA: binomial, n>p case", {

  library(SLOPE)
  set.seed(116)

  n = 100
  p = 10

  d <- randomProblem(n, p, response="binomial", density = 0.5)
  
  fista_solvers <- FISTA(d$x, d$y, family="binomial",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="binomial",alpha=c(1.0,0.005))
  
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: binomial, n<p case", {

  library(SLOPE)
  set.seed(114)
  n = 10
  p = 20

  d <- randomProblem(n, p, response="binomial", density = 0.5)

  fista_solvers <- FISTA(d$x, d$y, family="binomial",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="binomial",alpha=c(1.0,0.005))

  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: poisson, n>p case", {

  library(SLOPE)
  set.seed(164)

  n = 100
  p = 10

  d <- randomProblem(n, p, response="poisson", density = 0.5)
  
  fista_solvers <- FISTA(d$x, d$y, family="poisson",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="poisson",alpha=c(1.0,0.005))

  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: poisson, n<p case", {

  library(SLOPE)
  set.seed(24)

  n = 10
  p = 20

  d <- randomProblem(n, p, response="poisson", density = 0.5)
  
  fista_solvers <- FISTA(d$x, d$y, family="poisson",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="poisson",alpha=c(1.0,0.005))

  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: multinomial, n>p case", {
  library(SLOPE)
  set.seed(242)

  n = 100
  p = 10

  d <- randomProblem(n, p, response="multinomial", density = 0.5)
  
  fista_solvers <- FISTA(d$x, d$y, family="multinomial",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="multinomial",alpha=c(1.0,0.005))
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: multinomial, n<p case", {
  library(SLOPE)
  set.seed(514)

  n = 10
  p = 20

  d <- randomProblem(n, p, response="multinomial", density = 0.5)
  
  fista_solvers <- FISTA(d$x, d$y, family="multinomial",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="multinomial",alpha=c(1.0,0.005))
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})
