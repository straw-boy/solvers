test_that("Proximal Newton works same as current SLOPE package implementation", {
  # skip('Takes too much time')
  library(SLOPE)

  # fista_slope <- SLOPE(bodyfat$x, bodyfat$y,family="gaussian", solver = "fista")
  # pn_solvers <- PN(bodyfat$x, bodyfat$y,family="gaussian")
  # expect_equivalent(coef(fista_slope), coef(pn_solvers), tol = 1e-2)

  # fista_slope <- SLOPE(heart$x, heart$y,family="binomial", solver = "fista")
  # pn_solvers <- PN(heart$x, heart$y,family="binomial")
  # expect_equivalent(coef(fista_slope), coef(pn_solvers), tol = 1e-2)

  fista_slope <- ADMM(abalone$x, abalone$y,family="poisson", opt_algo = "nr")
  pn_solvers <- PN(abalone$x, abalone$y,family="poisson")
  expect_equivalent(coef(fista_slope), coef(pn_solvers), tol = 1e-2)
})


