# test_that("ADMM (as well as Newton-Raphson) works same as SLOPE package (excluding intercept) ", {
#   library(SLOPE)

#   fista_slope <- SLOPE(bodyfat$x, bodyfat$y, solver = "fista",intercept=FALSE)
#   admm_solvers <- ADMM(bodyfat$x, bodyfat$y,intercept=FALSE)
#   expect_equivalent(coef(fista_slope), coef(admm_solvers), tol = 1e-2)

#   fista_slope <- SLOPE(heart$x, heart$y,family="binomial", solver = "fista",intercept=FALSE)
#   admm_solvers <- ADMM(heart$x, heart$y,family="binomial",intercept=FALSE)
#   expect_equivalent(coef(fista_slope), coef(admm_solvers), tol = 1e-2)

#   fista_slope <- SLOPE(abalone$x, abalone$y,family="poisson", solver = "fista",intercept=FALSE)
#   admm_solvers <- ADMM(abalone$x, abalone$y,family="poisson",intercept=FALSE)
#   expect_equivalent(coef(fista_slope), coef(admm_solvers), tol = 1e-2)
# })

# test_that("ADMM (as well as Newton-Raphson) works (including intercept) ", {
#   library(SLOPE)

#   fista_solvers <- FISTA(bodyfat$x, bodyfat$y, solver = "fista")
#   admm_solvers <- ADMM(bodyfat$x, bodyfat$y)
#   expect_equivalent(coef(fista_solvers), coef(admm_solvers), tol = 1e-2)

#   fista_solvers <- FISTA(heart$x, heart$y,family="binomial", solver = "fista")
#   admm_solvers <- ADMM(heart$x, heart$y,family="binomial")
#   expect_equivalent(coef(fista_solvers), coef(admm_solvers), tol = 1e-2)

#   fista_solvers <- FISTA(abalone$x, abalone$y,family="poisson", solver = "fista")
#   admm_solvers <- ADMM(abalone$x, abalone$y,family="poisson")
#   expect_equivalent(coef(fista_solvers), coef(admm_solvers), tol = 1e-2)
# })


# test_that("ADMM (as well as Newton-Raphson) works same as SLOPE package (including intercept) ", {
#   skip('Intercept mismatch in this')

#   library(SLOPE)

#   fista_slope <- SLOPE(bodyfat$x, bodyfat$y, solver = "fista")
#   admm_solvers <- ADMM(bodyfat$x, bodyfat$y)
#   expect_equivalent(coef(fista_slope), coef(admm_solvers), tol = 1e-2)

#   fista_slope <- SLOPE(heart$x, heart$y,family="binomial", solver = "fista")
#   admm_solvers <- ADMM(heart$x, heart$y,family="binomial")
#   expect_equivalent(coef(fista_slope), coef(admm_solvers), tol = 1e-2)

#   fista_slope <- SLOPE(abalone$x, abalone$y,family="poisson", solver = "fista")
#   admm_solvers <- ADMM(abalone$x, abalone$y,family="poisson")
#   expect_equivalent(coef(fista_slope), coef(admm_solvers), tol = 1e-2)
# })

# test_that("ADMM output is same for different optimization algorithm choices (BFGS) ", {

#   admm_nr <- ADMM(bodyfat$x, bodyfat$y, opt_algo="nr")
#   admm_bfgs <- ADMM(bodyfat$x, bodyfat$y, opt_algo="bfgs")
#   expect_equivalent(coef(admm_nr), coef(admm_bfgs), tol = 1e-2)

#   admm_nr <- ADMM(heart$x, heart$y,family="binomial",opt_algo="nr")
#   admm_bfgs <- ADMM(heart$x, heart$y,family="binomial",opt_algo="bfgs")
#   expect_equivalent(coef(admm_nr), coef(admm_bfgs), tol = 1e-2)

#   admm_nr <- ADMM(abalone$x, abalone$y,family="poisson",opt_algo="nr")
#   admm_bfgs <- ADMM(abalone$x, abalone$y,family="poisson",opt_algo="bfgs")
#   expect_equivalent(coef(admm_nr), coef(admm_bfgs), tol = 1e-2)
# })

# test_that("ADMM output is same for different optimization algorithm choices (L-BFGS)", {

#   admm_nr <- ADMM(bodyfat$x, bodyfat$y, opt_algo="nr")
#   admm_lbfgs <- ADMM(bodyfat$x, bodyfat$y, opt_algo="lbfgs")
#   expect_equivalent(coef(admm_nr), coef(admm_lbfgs), tol = 1e-2)

#   admm_nr <- ADMM(heart$x, heart$y,family="binomial",opt_algo="nr")
#   admm_lbfgs <- ADMM(heart$x, heart$y,family="binomial",opt_algo="lbfgs")
#   expect_equivalent(coef(admm_nr), coef(admm_lbfgs), tol = 1e-2)

#   admm_nr <- ADMM(abalone$x, abalone$y,family="poisson",opt_algo="nr")
#   admm_lbfgs <- ADMM(abalone$x, abalone$y,family="poisson",opt_algo="lbfgs")
#   expect_equivalent(coef(admm_nr), coef(admm_lbfgs), tol = 1e-2)
# })

