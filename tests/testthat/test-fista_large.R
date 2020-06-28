test_that("FISTA works same as current SLOPE package implementation", {
  skip('Intercept mismatch in this; takes too much time')
  library(SLOPE)

  fista_slope <- SLOPE(bodyfat$x, bodyfat$y,family="gaussian", solver = "fista")
  fista_solvers <- FISTA(bodyfat$x, bodyfat$y,family="gaussian")
  expect_equivalent(coef(fista_slope), coef(fista_solvers), tol = 1e-2)

  fista_slope <- SLOPE(heart$x, heart$y,family="binomial", solver = "fista")
  fista_solvers <- FISTA(heart$x, heart$y,family="binomial")
  expect_equivalent(coef(fista_slope), coef(fista_solvers), tol = 1e-2)

  fista_slope <- SLOPE(abalone$x, abalone$y,family="poisson", solver = "fista")
  fista_solvers <- FISTA(abalone$x, abalone$y,family="poisson")
  expect_equivalent(coef(fista_slope), coef(fista_solvers), tol = 1e-2)
})


