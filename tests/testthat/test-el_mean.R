test_that("convergence check", {
  data("women")
  x <- women$weight
  grid <- seq(120, 160, length.out = 1000)
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  conv <- function(par) {
    el_mean(x, par, control = optcfg)@optim$convergence
  }
  expect_true(all(vapply(grid, conv, FUN.VALUE = logical(1))))
})

test_that("probabilities add up to 1", {
  data("women")
  x <- women$height
  par <- 65
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(x, par, control = optcfg)
  expect_output(print(fit))
  expect_equal(sum(exp(fit@logp)), 1, tolerance = 1e-07)
})

test_that("probabilities add up to 1 (weighted)", {
  data("women")
  x <- women$height
  par <- 65
  w <- women$weight
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(x, par, weights = w, optcfg)
  expect_equal(sum(exp(fit@logp)), 1, tolerance = 1e-7)
})

test_that("identical weights == no weights", {
  data("women")
  x <- women$height
  par <- 65
  w <- rep(0.5, length(x))
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit1 <- el_mean(x, par, control = optcfg)
  fit2 <- el_mean(x, par, weights = w, control = optcfg)
  fit2@weights <- numeric()
  expect_equal(fit1, fit2)
})

test_that("loglik to loglr", {
  data("sleep")
  x <- sleep$extra
  n <- length(x)
  par <- 0
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(x, par, control = optcfg)
  expect_equal(fit@logl + n * log(n), fit@loglr)
})

test_that("loglik to loglr (weighted)", {
  data("women")
  x <- women$height
  n <- length(x)
  par <- 65
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(x, par, weights = women$weight, control = optcfg)
  w <- weights(fit)
  expect_equal(fit@logl + sum(w * (log(n) - log(w))), fit@loglr)
})

test_that("non-full rank", {
  x <- matrix(c(1, 1, 2, 2), ncol = 2)
  w <- c(1, 2)
  par <- c(0, 0)
  optcfg <- el_control(maxit_l = 20L, tol_l = 1e-08, th = 1e+10)
  expect_error(el_mean(par, x, control = optcfg))
  expect_error(el_mean(par, x, weights = w, control = optcfg))
})

test_that("invalid 'x", {
  par <- 0
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  expect_error(el_mean(c(1, Inf), par, control = optcfg))
  expect_error(el_mean(10, par, control = optcfg))
})

test_that("invalid 'par", {
  data("sleep")
  x <- sleep$extra
  par <- 0
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  expect_error(el_mean(x, NA, control = optcfg))
  expect_error(el_mean(x, NULL, control = optcfg))
  expect_error(el_mean(x, Inf, control = optcfg))
  expect_error(el_mean(x, NaN, control = optcfg))
  expect_error(el_mean(x, c(0, 0), control = optcfg))
})

test_that("identical results for repeated executions", {
  data("women")
  x <- women$height
  par <- 65
  w <- women$weight
  fit1 <- el_mean(x, par, control = el_control(th = 1e+10))
  fit2 <- el_mean(x, par, control = el_control(th = 1e+10))
  expect_equal(fit1, fit2)

  wfit1 <- el_mean(x, par, weights = w, control = el_control(th = 1e+10))
  wfit2 <- el_mean(x, par, weights = w, control = el_control(th = 1e+10))
  expect_equal(wfit1, wfit2)

  lhs <- 4
  elt1 <- elt(fit1, lhs = lhs)
  elt2 <- elt(fit2, lhs = lhs)
  expect_equal(elt1, elt2)

  elt3 <- elt(wfit1, lhs = lhs)
  elt4 <- elt(wfit2, lhs = lhs)
  expect_equal(elt3, elt4)
})
