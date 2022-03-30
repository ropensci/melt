test_that("invalid 'level'", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- (max(x) + min(x)) / 2
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  expect_error(confint(fit, level = -1))
  expect_error(confint(fit, level = 2))
  expect_error(confint(fit, level = Inf))
  expect_error(confint(fit, level = c(0, 0)))
})

test_that("invalid 'parm'", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- (max(x) + min(x)) / 2
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  expect_error(confint(fit, parm = NA))
  expect_error(confint(fit, parm = NULL))
  expect_error(confint(fit, parm = NaN))
  expect_error(confint(fit, parm = Inf))
})

test_that("'level' == 1", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- (max(x) + min(x)) / 2
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  ci <- confint(fit, level = 1)
  expect_equal(ci[1], -Inf)
  expect_equal(ci[2], Inf)
})

test_that("'level' == 0", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- (max(x) + min(x)) / 2
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  ci <- confint(fit, level = 0)
  expect_equal(ci[2] - ci[1], 0)
})

test_that("empty model", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  y <- 1 + x + rnorm(n)
  df <- data.frame(y, x)
  fit <- el_lm(y~ -1, df)
  ci <- confint(fit)
  expect_equal(nrow(ci), 0)
})
