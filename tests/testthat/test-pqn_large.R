test_that("Proximal Newton works same as ADMM", {
  skip('Takes too much time')
  library(SLOPE)

  pn_solvers <- PN(bodyfat$x, bodyfat$y, family="gaussian",hessian_calc="exact")
  pqn_solvers <- PN(bodyfat$x, bodyfat$y, family="gaussian",hessian_calc="lbfgs")
  
  expect_equivalent(coef(pn_solvers),coef(pqn_solvers), tol = 1e-2)

  pn_solvers <- PN(heart$x, heart$y, family="binomial",hessian_calc="exact")
  pqn_solvers <- PN(heart$x, heart$y, family="binomial",hessian_calc="lbfgs")
  expect_equivalent(coef(pn_solvers),coef(pqn_solvers), tol = 1e-2)

  pn_solvers <- PN(abalone$x, abalone$y, family="poisson",hessian_calc="exact")
  pqn_solvers <- PN(abalone$x, abalone$y, family="poisson",hessian_calc="lbfgs")
  expect_equivalent(coef(pn_solvers),coef(pqn_solvers), tol = 1e-2)
})


