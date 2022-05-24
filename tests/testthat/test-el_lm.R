test_that("probabilities add up to 1", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  # expect_output(print(formula(fit)))
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_equal(sum(exp(fit@logp)), 1)
})

test_that("probabilities add up to 1 (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, weights = w, control = optcfg)
  expect_equal(sum(exp(fit@logp)), 1)
})

test_that("loglik to loglr", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  expect_equal(fit@logl + n * log(n), fit@loglr)
})

test_that("loglik to loglr (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, weights = w, control = optcfg)
  w <- fit@weights
  expect_equal(fit@logl + sum(w * (log(n) - log(w))), fit@loglr)
})

test_that("non-full rank", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- x
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  expect_error(el_lm(y ~ x + x2, df, control = optcfg))
  expect_error(el_lm(y ~ x + x2, df, weights = w, control = optcfg))
})

test_that("empty model", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  y <- 1 + x + rnorm(n)
  df <- data.frame(y, x)
  fit <- el_lm(y ~ 0, df)
  expect_output(print(summary(fit)))
})

test_that("same results with parallel computing", {
  skip_on_os("windows", arch = "i386")
  n <- 500
  p <- 15
  b <- rnorm(p)
  x <- matrix(rnorm(n * p), ncol = p)
  y <- 1 + x %*% as.vector(b) + rnorm(n)
  df <- data.frame(y, x)
  fit <- el_lm(y ~ ., df, control = el_control(th = 1e+10, nthreads = 1))
  fit2 <- el_lm(y ~ ., df, control = el_control(th = 1e+10))
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)

  w <- 1 + runif(n, min = -0.5, max = 0.5)
  wfit <- el_lm(y ~ ., df, weights = w,
               control = el_control(th = 1e+10, nthreads = 1))
  wfit2 <- el_lm(y ~ ., df, weights = w,
                control = el_control(th = 1e+10))
  expect_equal(wfit@optim, wfit2@optim)
  expect_equal(wfit@parTests, wfit2@parTests)
})
