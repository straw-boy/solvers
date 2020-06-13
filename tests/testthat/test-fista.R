
test_that("FISTA works same as current SLOPE package implementation", {
  skip("aa");
  fista_slope <- SLOPE(bodyfat$x, bodyfat$y, solver = "fista")
  fista_solvers <- FISTA(bodyfat$x, bodyfat$y)
  print(coef(fista_slope)-coef(fista_solvers))
  expect_equivalent(coef(fista_slope), coef(fista_solvers), tol = 1e-3)
})
