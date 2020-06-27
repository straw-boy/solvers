test_that("FISTA: gaussian, n>p case", {
  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- randomProblem(n,p,response="gaussian")
  
  fista_solvers <- FISTA(d$x, d$y, family="gaussian",path_length=3)
  fista_slope <- SLOPE(d$x, d$y, family="gaussian",path_length=3)
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: gaussian, n<p case", {
  library(SLOPE)
  set.seed(1)


  n = 10
  p = 20

  d <- randomProblem(n,p,response="gaussian")
  
  fista_solvers <- FISTA(d$x, d$y, family="gaussian",path_length=3)
  fista_slope <- SLOPE(d$x, d$y, family="gaussian",path_length=3)
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})


test_that("FISTA: binomial, n>p case", {
  # skip('Intercept mismatch in this')

  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- randomProblem(n,p,response="binomial")
  
  fista_solvers <- FISTA(d$x, d$y, family="binomial",path_length=3)
  fista_slope <- SLOPE(d$x, d$y, family="binomial",path_length=3)
  
  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

test_that("FISTA: binomial, n<p case", {
  skip('Intercept mismatch in this')

  library(SLOPE)
  set.seed(1)
  n = 10
  p = 20

  d <- randomProblem(n,p,response="binomial")
  
  fista_solvers <- FISTA(d$x, d$y, family="binomial",path_length=3)
  fista_slope <- SLOPE(d$x, d$y, family="binomial",path_length=3)

  expect_equivalent(coef(fista_solvers), coef(fista_slope), tol = 1e-2)

})

