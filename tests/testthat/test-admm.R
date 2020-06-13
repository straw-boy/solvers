test_that("ADMM (as well as Newton-Raphson) works", {
  fista_slope <- SLOPE(bodyfat$x, bodyfat$y, solver = "fista")
  admm_solvers <- ADMM(bodyfat$x, bodyfat$y)

  expect_equivalent(coef(fista_slope), coef(admm_solvers), tol = 1e-3)
})
