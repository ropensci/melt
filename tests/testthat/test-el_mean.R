test_that("convergence check", {
  skip_on_os("windows", arch = "i386")
  x <- c(-1.5, 1.5, rnorm(10))
  grid <- seq(-1, 1, length.out = 1000)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  conv <- function(par) {
    el_mean(par, x, control = optcfg)$optim$convergence
  }
  expect_true(all(vapply(grid, conv, FUN.VALUE = logical(1))))
})

test_that("probabilities add up to 1", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- (max(x) + min(x)) / 2
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  expect_output(print(fit))
  expect_equal(sum(exp(fit$log.prob)), 1)
})

test_that("probabilities add up to 1 (weighted)", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- (max(x) + min(x)) / 2
  w <- 1 + runif(10, min = -0.5, max = 0.5)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, w, optcfg)
  expect_equal(sum(exp(fit$log.prob)), 1)
})

test_that("identical weights == no weights", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  w <- rep(runif(1), length(x))
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  a1 <- el_mean(par, x, control = optcfg)$optim
  a2 <- el_mean(par, x, w, control = optcfg)$optim
  expect_equal(a1$lambda, a2$lambda)
  expect_equal(a1$logLR, a2$logLR)
})

test_that("loglik to loglr", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- (max(x) + min(x)) / 2
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  expect_equal(fit$loglik + n * log(n), fit$optim$logLR)
})

test_that("loglik to loglr (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- (max(x) + min(x)) / 2
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, w, control = optcfg)
  w <- fit$weights
  expect_equal(fit$loglik + sum(w * (log(n) - log(w))), fit$optim$logLR)
})

test_that("non-full rank", {
  skip_on_os("windows", arch = "i386")
  x <- matrix(c(1, 1, 2, 2), ncol = 2)
  w <- c(1, 2)
  par <- c(0, 0)
  optcfg <- list(maxit = 20L, tol = 1e-08, th = 1e+10)
  expect_error(el_mean(par, x, control = optcfg))
  expect_error(el_mean(par, x, w, control = optcfg))
})

test_that("invalid 'x", {
  skip_on_os("windows", arch = "i386")
  par <- 0
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  expect_error(el_mean(par, c(1, Inf), control = optcfg))
  expect_error(el_mean(par, rnorm(1), control = optcfg))
})

test_that("invalid 'par", {
  skip_on_os("windows", arch = "i386")
  x <- matrix(c(1, 1, 2, 2), ncol = 2)
  par <- 0
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  expect_error(el_mean(par, x, control = optcfg))
  expect_error(el_mean(NA, rnorm(10), control = optcfg))
  expect_error(el_mean(NULL, rnorm(10), control = optcfg))
  expect_error(el_mean(Inf, rnorm(10), control = optcfg))
  expect_error(el_mean(NaN, rnorm(10), control = optcfg))
})
