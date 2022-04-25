test_that("convergence check", {
  skip_on_os("windows", arch = "i386")
  x <- c(-1.5, 1.5, rnorm(10))
  grid <- seq(-1, 1, length.out = 1000)
  optcfg <- control_el(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  conv <- function(par) {
    el_mean2(par, x, control = optcfg)@optim$convergence
  }
  expect_true(all(vapply(grid, conv, FUN.VALUE = logical(1))))
})

test_that("probabilities add up to 1", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  optcfg <- control_el(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  expect_output(print(fit))
  expect_equal(sum(exp(fit$log.prob)), 1)
})

test_that("probabilities add up to 1 (weighted)", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  w <- 1 + runif(10, min = -0.5, max = 0.5)
  optcfg <- control_el(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean2(par, x, weights = w, optcfg)
  expect_equal(sum(exp(fit@logp)), 1)
})

test_that("identical weights == no weights", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  w <- rep(runif(1), length(x))
  optcfg <- control_el(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  a1 <- el_mean2(par, x, control = optcfg)
  a2 <- el_mean2(par, x, weights = w, control = optcfg)
  a2@weights <- w
  expect_equal(a1, a2)
})

test_that("loglik to loglr", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- runif(1, min(x), max(x))
  optcfg <- control_el(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean2(par, x, control = optcfg)
  expect_equal(fit@logl + n * log(n), fit@optim$logLR)
})

test_that("loglik to loglr (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- runif(1, min(x), max(x))
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- control_el(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean2(par, x, weights = w, control = optcfg)
  w <- fit@weights
  expect_equal(fit@logl + sum(w * (log(n) - log(w))), fit@optim$logLR)
})

test_that("non-full rank", {
  skip_on_os("windows", arch = "i386")
  x <- matrix(c(1, 1, 2, 2), ncol = 2)
  w <- c(1, 2)
  par <- c(0, 0)
  optcfg <- control_el(maxit_l = 20L, tol_l = 1e-08, th = 1e+10)
  expect_error(el_mean2(par, x, control = optcfg))
  expect_error(el_mean2(par, x, weights = w, control = optcfg))
})

test_that("invalid 'x", {
  skip_on_os("windows", arch = "i386")
  par <- 0
  optcfg <- control_el(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  expect_error(el_mean2(par, c(1, Inf), control = optcfg))
  expect_error(el_mean2(par, rnorm(1), control = optcfg))
})

test_that("invalid 'par", {
  skip_on_os("windows", arch = "i386")
  x <- matrix(c(1, 1, 2, 2), ncol = 2L)
  par <- 0
  optcfg <- control_el(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  expect_error(el_mean2(par, x, control = optcfg))
  expect_error(el_mean2(NA, rnorm(10), control = optcfg))
  expect_error(el_mean2(NULL, rnorm(10), control = optcfg))
  expect_error(el_mean2(Inf, rnorm(10), control = optcfg))
  expect_error(el_mean2(NaN, rnorm(10), control = optcfg))
  x <- rnorm(10)
  par <- c(0, 0)
  expect_error(el_mean2(par, x, control = optcfg))
  expect_error(el_mean2(0, x, control = list(maxit = 200)))
})
