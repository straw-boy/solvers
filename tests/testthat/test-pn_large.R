test_that("Proximal Newton works same as ADMM", {
  skip('Takes too much time')
  library(SLOPE)

  admm_solvers <- ADMM(bodyfat$x, bodyfat$y, family="gaussian", opt_algo = "nr")
  pn_solvers <- PN(bodyfat$x, bodyfat$y, family="gaussian")
  
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

  admm_solvers <- ADMM(heart$x, heart$y, family="binomial", opt_algo="nr")
  pn_solvers <- PN(heart$x, heart$y, family="binomial")
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

  admm_solvers <- ADMM(abalone$x, abalone$y, family="poisson", opt_algo = "nr")
  pn_solvers <- PN(abalone$x, abalone$y, family="poisson")
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)

  admm_solvers <- ADMM(wine$x, wine$y, family="multinomial", opt_algo = "nr")
  pn_solvers <- PN(wine$x, wine$y, family="multinomial")
  expect_equivalent(coef(admm_solvers), coef(pn_solvers), tol = 1e-2)
})


