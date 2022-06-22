test_that("invalid 'object", {
  data("mtcars")
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_eval(mtcars$mpg, control = optcfg)
  expect_error(elt(fit, lhs = 1, control = optcfg))
  fit2 <- el_lm(mpg ~ 0, mtcars)
  expect_error(elt(fit2, lhs = 1, control = optcfg))
})

test_that("invalid 'lhs' and 'rhs'", {
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x + x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  # invalid rhs
  lhs <- matrix(c(0, 1, -1), nrow = 1)
  expect_error(elt(fit, control = optcfg))
  expect_error(elt(fit, rhs = c(NA, 0), lhs = lhs, control = optcfg))
  expect_error(elt(fit, rhs = c(1, 0), lhs = lhs, control = optcfg))
  expect_error(elt(fit, rhs = matrix(c(1, 0, 0), ncol = 3), lhs = lhs,
                   control = optcfg))
  expect_error(elt(fit, rhs = matrix(c(1, 0, 0), ncol = 1), lhs = lhs,
                   control = optcfg))
  expect_error(elt(fit, rhs = matrix(c("error"), ncol = 1), lhs = lhs,
                   control = optcfg))

  w <- 1 + runif(n, min = -0.5, max = 0.5)
  fit2 <- el_lm(y ~ x + x2, df, weights = w, control = optcfg)
  out <- elt(fit2, lhs = lhs, control = optcfg)
  expect_lt(out@statistic, 20)

  # invalid lhs
  lhs <- matrix(c(1, 0, 0, 0, NA, 1), nrow = 2)
  expect_error(elt(fit, lhs = lhs, control = optcfg))
  expect_error(elt(fit2, lhs = lhs, control = optcfg))
  lhs <- matrix(rnorm(4), ncol = 2)
  expect_error(elt(fit, lhs = lhs, control = optcfg))
  expect_error(elt(fit2, lhs = lhs, control = optcfg))
  lhs <- matrix(rnorm(12), ncol = 3)
  expect_error(elt(fit, lhs = lhs, control = optcfg))
  expect_error(elt(fit2, lhs = lhs, control = optcfg))
  lhs <- matrix(c(1, 1, 0, 0, 0, 0), nrow = 2)
  expect_error(elt(fit, lhs = lhs, control = optcfg))
  expect_error(elt(fit2, lhs = lhs, control = optcfg))
})

test_that("when elt == eval", {
  skip_on_os("windows", arch = "i386")
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  # mean
  n <- 100
  x <- rnorm(n)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  par <- runif(1, min(x), max(x))
  fit <- el_mean(x, par, control = optcfg)
  fit2 <- elt(fit, rhs = par, control = optcfg)
  expect_equal(fit@optim$lambda, fit2@optim$lambda)
  # mean (weighted)
  fit <- el_mean(x, par, weights = w, control = optcfg)
  fit2 <- elt(fit, rhs = par, control = optcfg)
  expect_equal(fit@optim$lambda, fit2@optim$lambda)
  # lm
  x2 <- rnorm(n)
  y <- 1 + x + x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  lhs <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 2)
  rhs <- c(0, 0)
  fit2 <- elt(fit, lhs = lhs, rhs = rhs, control = optcfg)
  expect_output(print(fit2))
  expect_equal(fit@optim$lambda, fit2@optim$lambda)
  # lm (weighted)
  fit <- el_lm(y ~ x + x2, df, weights = w, control = optcfg)
  lhs <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 2)
  rhs <- c(0, 0)
  fit2 <- elt(fit, lhs = lhs, rhs = rhs, control = optcfg)
  expect_equal(fit@optim$lambda, fit2@optim$lambda)
})

test_that("invalid 'calibrate'", {
  data("mtcars")
  fit <- el_lm(mpg ~ wt, data = mtcars)
  expect_error(elt(fit, rhs = c(1, 1), calibrate = "f"))
})

test_that("invalid 'control'", {
  data("sleep")
  fit <- el_mean(sleep$extra, 0)
  expect_error(elt(fit, lhs = 1, control = list(maxit = 200L)))
})

test_that("vector 'lhs'", {
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x + x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  expect_error(elt(fit, lhs = c(1, -1, 0, 0)))
  out <- elt(fit, lhs = c(1, -1, 0))
  expect_s4_class(out, "ELT")
  expect_output(show(out))
  expect_output(print(out))
})

test_that("calibration", {
  data("women")
  fit <- el_mean(women$height, 65)
  out <- elt(fit, rhs = 67, calibrate = "f")
  out2 <- elt(fit, rhs = 67, calibrate = "boot")
  expect_s4_class(out, "ELT")
  expect_s4_class(out2, "ELT")
})
