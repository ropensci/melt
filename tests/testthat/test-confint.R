test_that("invalid 'level'", {
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  expect_error(confint(fit, level = -1))
  expect_error(confint(fit, level = 2))
  expect_error(confint(fit, level = Inf))
  expect_error(confint(fit, level = c(0, 0)))
})

test_that("invalid 'parm'", {
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  expect_error(confint(fit, parm = NA))
  expect_error(confint(fit, parm = NULL))
  expect_error(confint(fit, parm = NaN))
  expect_error(confint(fit, parm = Inf))
})

test_that("'level' == 1", {
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  ci <- confint(fit, level = 1)
  expect_equal(ci[1], -Inf)
  expect_equal(ci[2], Inf)
})

test_that("'level' == 0", {
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  ci <- confint(fit, level = 0)
  expect_equal(ci[2] - ci[1], 0)
})

test_that("empty model", {
  n <- 100
  x <- rnorm(n)
  y <- 1 + x + rnorm(n)
  df <- data.frame(y, x)
  fit <- el_lm(y ~ -1, df)
  ci <- confint(fit)
  expect_equal(nrow(ci), 0)
})

test_that("no effect of nthreads", {
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 2 * x + 3 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  fit <- el_lm(y ~ x + x2)
  parm <- sample(c(1, 2, 3), size = 10, replace = TRUE)
  ci1 <- confint(fit, parm = parm, control = el_control(nthreads = 1L))
  ci2 <- confint(fit, parm = parm)
  expect_equal(ci1, ci2)
})
