test_that("FISTA: gaussian, n>p case", {
  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="gaussian")
  
  fista_solvers <- FISTA(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: gaussian, n<p case", {
  library(SLOPE)
  set.seed(1)


  n = 5
  p = 20

  d <- randomProblem(n,p,response="gaussian")
  
  fista_solvers <- FISTA(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="gaussian",alpha=c(1.0,0.005))
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})


test_that("FISTA: binomial, n>p case", {
  # skip('Intercept mismatch in this')

  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="binomial")
  
  fista_solvers <- FISTA(d$x, d$y, family="binomial",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="binomial",alpha=c(1.0,0.005))
  
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: binomial, n<p case", {
  skip('Intercept mismatch in this')

  library(SLOPE)
  set.seed(1)
  n = 10
  p = 20

  d <- randomProblem(n,p,response="binomial")
  print(d)
  
  fista_solvers <- FISTA(d$x, d$y, family="binomial",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="binomial",alpha=c(1.0,0.005))

  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: poisson, n>p case", {
  skip('Intercept mismatch in this')

  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="poisson")
  
  fista_solvers <- FISTA(d$x, d$y, family="poisson",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="poisson",alpha=c(1.0,0.005))

  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: poisson, n<p case", {
  skip('Intercept mismatch in this')

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="poisson")
  
  fista_solvers <- FISTA(d$x, d$y, family="poisson",alpha=c(1.0,0.005))
  fista_slope <- SLOPE(d$x, d$y, family="poisson",alpha=c(1.0,0.005))

  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})
