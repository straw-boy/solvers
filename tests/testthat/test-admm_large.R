
test_that("ADMM works same as FISTA", {
  
  # skip("Takes too much time")
  library(SLOPE)

  fista_solvers <- FISTA(bodyfat$x, bodyfat$y)
  admm_solvers <- ADMM(bodyfat$x, bodyfat$y, opt_algo = "nr")
  expect_equivalent(coef(fista_solvers), coef(admm_solvers), tol = 1e-2)

  fista_solvers <- FISTA(heart$x, heart$y,family="binomial")
  admm_solvers <- ADMM(heart$x, heart$y,family="binomial", opt_algo = "nr")
  expect_equivalent(coef(fista_solvers), coef(admm_solvers), tol = 1e-2)

  fista_solvers <- FISTA(abalone$x, abalone$y,family="poisson")
  admm_solvers <- ADMM(abalone$x, abalone$y,family="poisson", opt_algo = "nr")
  expect_equivalent(coef(fista_solvers), coef(admm_solvers), tol = 1e-2)

  fista_solvers <- FISTA(wine$x, wine$y,family="multinomial")
  admm_solvers <- ADMM(wine$x, wine$y,family="multinomial",opt_algo="nr")
  expect_equivalent(coef(fista_solvers), coef(admm_solvers), tol = 1e-2)

})

test_that("ADMM output is same for different optimization algorithm choices", {

  skip("Takes too much time")
  library(SLOPE)

  admm_nr <- ADMM(bodyfat$x, bodyfat$y, family="gaussian", opt_algo="nr")
  admm_bfgs <- ADMM(bodyfat$x, bodyfat$y, family="gaussian", opt_algo="bfgs")
  expect_equivalent(coef(admm_nr), coef(admm_bfgs), tol = 1e-2)

  admm_lbfgs <- ADMM(bodyfat$x, bodyfat$y, opt_algo="lbfgs")
  expect_equivalent(coef(admm_nr), coef(admm_lbfgs), tol = 1e-2)

  admm_nr <- ADMM(heart$x, heart$y, family="binomial",opt_algo="nr")
  admm_bfgs <- ADMM(heart$x, heart$y, family="binomial",opt_algo="bfgs")
  expect_equivalent(coef(admm_nr), coef(admm_bfgs), tol = 1e-2)

  admm_lbfgs <- ADMM(heart$x, heart$y, family="binomial",opt_algo="lbfgs")
  expect_equivalent(coef(admm_nr), coef(admm_lbfgs), tol = 1e-2)

  admm_nr <- ADMM(abalone$x, abalone$y, family="poisson",opt_algo="nr")
  admm_bfgs <- ADMM(abalone$x, abalone$y, family="poisson",opt_algo="bfgs")
  expect_equivalent(coef(admm_nr), coef(admm_bfgs), tol = 1e-2)

  admm_lbfgs <- ADMM(abalone$x, abalone$y, family="poisson",opt_algo="lbfgs")
  expect_equivalent(coef(admm_nr), coef(admm_lbfgs), tol = 1e-2)

  admm_nr <- ADMM(wine$x, wine$y, family="multinomial",opt_algo="nr")
  admm_bfgs <- ADMM(wine$x, wine$y, family="multinomial",opt_algo="bfgs")
  expect_equivalent(coef(admm_nr), coef(admm_bfgs), tol = 1e-2)

  admm_lbfgs <- ADMM(wine$x, wine$y, family="multinomial",opt_algo="lbfgs")
  expect_equivalent(coef(admm_nr), coef(admm_lbfgs), tol = 1e-2)

})
