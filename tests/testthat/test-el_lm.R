test_that("probabilities add up to 1", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  expect_output(print(formula(fit)))
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_equal(sum(exp(fit$log.prob)), 1)
})

test_that("probabilities add up to 1 (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, weights =  w, control = optcfg)
  expect_equal(sum(exp(fit$log.prob)), 1)
})

test_that("loglik to loglr", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  expect_equal(fit$loglik + n * log(n), fit$optim$logLR)
})

test_that("loglik to loglr (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, weights =  w, control = optcfg)
  w <- fit$weights
  expect_equal(fit$loglik + sum(w * (log(n) - log(w))), fit$optim$logLR)
})

test_that("non-full rank", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- x
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  expect_error(el_lm(y ~ x + x2, df, control = optcfg))
  expect_error(el_lm(y ~ x + x2, df, weights =  w, control = optcfg))
})
