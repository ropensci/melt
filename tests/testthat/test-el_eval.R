test_that("convergence check", {
  skip_on_os("windows", arch = "i386")
  x <- c(-1.5, 1.5, rnorm(10))
  grid <- seq(-1, 1, length.out = 1000)
  conv <- function(par) {
    el_eval(x - par,
      control =
        el_control(
          maxit_l = 20L, tol_l = 1e-08,
          th = 1e+10
        )
    )$optim$convergence
  }
  expect_true(all(vapply(grid, conv, FUN.VALUE = logical(1))))
})

test_that("probabilities add up to 1", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_eval(x - par, control = optcfg)
  expect_equal(sum(exp(fit$logp)), 1)
})

test_that("probabilities add up to 1 (weighted)", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  w <- 1 + runif(10, min = -0.5, max = 0.5)
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_eval(x - par, w, optcfg)
  expect_equal(sum(exp(fit$logp)), 1)
})

test_that("loglik to loglr", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_eval(x - par, control = optcfg)
  expect_equal(fit$logl + n * log(n), fit$loglr)
})

test_that("loglik to loglr (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- runif(1, min(x), max(x))
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_eval(x - par, w, optcfg)
  w <- fit$weights
  expect_equal(fit$logl + sum(w * (log(n) - log(w))), fit$loglr)
})

test_that("identical weights == no weights", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  g <- x - par
  w <- rep(runif(1), length(x))
  optcfg <- el_control(maxit_l = 20L, tol_l = 1e-08, th = 1e+10)
  a1 <- el_mean(par, x, control = optcfg)
  a2 <- el_mean(par, x, weights = w, control = optcfg)
  expect_equal(a1@optim, a2@optim)
})

test_that("non-full rank", {
  skip_on_os("windows", arch = "i386")
  g <- matrix(c(1, 1, 2, 2), ncol = 2)
  w <- c(1, 2)
  optcfg <- el_control(maxit_l = 20L, tol_l = 1e-08, th = 1e+10)
  expect_error(el_eval(g, control = optcfg))
  expect_error(el_eval(g, w, control = optcfg))
})

test_that("invalid 'g'", {
  skip_on_os("windows", arch = "i386")
  optcfg <- el_control(maxit_l = 20L, tol_l = 1e-08, th = 1e+10)
  expect_error(el_eval(matrix(rnorm(2), ncol = 2), control = optcfg))
  expect_error(el_eval(matrix(c(1, 1, 2, NA), ncol = 2), control = optcfg))
})

test_that("invalid 'control'", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  expect_error(el_eval(x - par, control = list(maxit = 200L)))
})
