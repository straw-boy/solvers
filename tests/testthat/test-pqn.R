test_that("Proximal Quasi-Newton: gaussian, n>p case", {
  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- solvers::randomProblem(n, p, response="gaussian", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), opt_algo="nr")
  pqn_solvers <- PN(d$x, d$y, family="gaussian", alpha=c(1.0,0.005), hessian_calc="lbfgs")
  expect_equivalent(coef(admm_solvers), coef(pqn_solvers), tol = 1e-2)

})

test_that("Proximal Quasi-Newton: gaussian, n<p case", {

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- solvers::randomProblem(n, p, response="gaussian", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="gaussian",alpha=c(1.0,0.005),opt_algo="nr")
  pqn_solvers <- PN(d$x, d$y, family="gaussian",alpha=c(1.0,0.005), hessian_calc="lbfgs")
  expect_equivalent(coef(admm_solvers), coef(pqn_solvers), tol = 1e-2)

})

test_that("Proximal Quasi-Newton: binomial, n>p case", {
  
  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- solvers::randomProblem(n, p, response="binomial", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005),opt_algo="nr")
  pqn_solvers <- PN(d$x, d$y, family="binomial",alpha=c(1.0,0.005), hessian_calc="lbfgs")
  expect_equivalent(coef(admm_solvers), coef(pqn_solvers), tol = 1e-2)

})


test_that("Proximal Quasi-Newton: binomial, n<p case", {
  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- solvers::randomProblem(n, p, response="binomial", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="binomial",alpha=c(1.0,0.005),opt_algo="nr", tol_abs=0,tol_rel=0,max_passes=10000)
  pqn_solvers <- PN(d$x, d$y, family="binomial",alpha=c(1.0,0.005), hessian_calc="lbfgs")
  expect_equivalent(coef(admm_solvers), coef(pqn_solvers), tol = 1e-2)

})

test_that("Proximal Quasi-Newton: poisson, n>p case", {

  library(SLOPE)
  set.seed(1)

  n = 100
  p = 10

  d <- solvers::randomProblem(n, p, response="poisson", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005),opt_algo="nr")
  pqn_solvers <- PN(d$x, d$y, family="poisson",alpha=c(1.0,0.005), hessian_calc="lbfgs")
  expect_equivalent(coef(admm_solvers), coef(pqn_solvers), tol = 1e-2)

})

test_that("Proximal Quasi-Newton: poisson, n<p case", {

  library(SLOPE)
  set.seed(1)

  n = 10
  p = 20

  d <- solvers::randomProblem(n, p, response="poisson", density = 0.5)
  
  admm_solvers <- ADMM(d$x, d$y, family="poisson",alpha=c(1.0,0.005),opt_algo="nr")
  # Passes if call is : PN(d$x, d$y, family="poisson",alpha=c(1.0,0.005), hessian_calc="lbfgs",tol=0, max_passes=550)
  pqn_solvers <- PN(d$x, d$y, family="poisson",alpha=c(1.0,0.005), hessian_calc="lbfgs")
  expect_equivalent(coef(admm_solvers), coef(pqn_solvers), tol = 1e-2)

})

