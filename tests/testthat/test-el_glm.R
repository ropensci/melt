test_that("probabilities add up to 1", {
  skip_on_os("windows", arch = "i386")
  n <- 50
  x <- rnorm(n)
  x2 <- rnorm(n)
  l <- -2 + 0.2 * x + 1 * x2
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_glm(y ~ x + x2, family = binomial, df, control = optcfg)
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_equal(sum(exp(fit@logp)), 1)
})

test_that("probabilities add up to 1 (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 50
  x <- rnorm(n)
  x2 <- rnorm(n)
  l <- -2 + 0.2 * x + 1 * x2
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  df <- data.frame(y, x, x2)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_glm(y ~ x + x2, family = binomial, df, control = optcfg)
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_equal(sum(exp(fit@logp)), 1)
})

test_that("loglik to loglr", {
  skip_on_os("windows", arch = "i386")
  n <- 50
  x <- rnorm(n)
  x2 <- rnorm(n)
  l <- -1 + 0.9 * x + 0.3 * x2
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_glm(y ~ x + x2, family = binomial, df, control = optcfg)
  expect_equal(fit@logl + n * log(n), fit@loglr)
})

test_that("loglik to loglr (weighted)", {
  skip_on_os("windows", arch = "i386")
  n <- 50
  x <- rnorm(n)
  x2 <- rnorm(n)
  l <- -1 + 0.9 * x + 0.3 * x2
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  df <- data.frame(y, x, x2)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- suppressWarnings(el_glm(y ~ x + x2, family = binomial, df,
                                  weights = w, control = optcfg))
  w <- fit@weights
  expect_equal(fit@logl + sum(w * (log(n) - log(w))), fit@loglr)
})

test_that("non-full rank", {
  skip_on_os("windows", arch = "i386")
  n <- 50
  x <- rnorm(n)
  x2 <- x
  l <- -1 + 0.9 * x + 0.3 * x2
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  expect_error(el_glm(y ~ x + x2, family = binomial, df, control = optcfg))
})

test_that("empty model", {
  skip_on_os("windows", arch = "i386")
  n <- 50
  x <- rnorm(n)
  x2 <- x
  l <- -1 + 0.9 * x + 0.3 * x2
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_glm(y ~ 0, family = binomial, control = optcfg)
  expect_output(print(summary(fit)))
})

test_that("invalid 'family'", {
  skip_on_os("windows", arch = "i386")
  n <- 50
  x <- rnorm(n)
  y <- -1 + 0.9 * x + rnorm(n)
  df <- data.frame(y, x)
  invalid <- function() {}
  expect_error(capture.output(el_glm(y ~ x, family = invalid, df)))
})

test_that("invalid 'control'", {
  skip_on_os("windows", arch = "i386")
  n <- 50
  x <- rnorm(n)
  x2 <- x
  l <- -1 + 0.9 * x + 0.3 * x2
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  df <- data.frame(y, x, x2)
  expect_error(el_glm(y ~ x + x2,
    family = binomial, df,
    control = list(maxit = 2L)
  ))
})

test_that("same results with parallel computing (binomial)", {
  skip_on_os("windows", arch = "i386")
  # skip_on_ci()
  n <- 500
  p <- 15
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  l <- 1 + x %*% as.vector(b)
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  df <- data.frame(cbind(y, x))
  lfit <- el_glm(y ~ ., family = binomial(link = "logit"), df,
                control = el_control(tol = 1e-08, th = 1e+10, nthreads = 1))
  lfit2 <- el_glm(y ~ ., family = binomial(link = "logit"), df,
                 control = el_control(tol = 1e-08, th = 1e+10))
  expect_equal(lfit@optim, lfit2@optim)
  expect_equal(lfit@parTests, lfit2@parTests)
  wlfit <- suppressWarnings(
    el_glm(y ~ ., family = binomial(link = "logit"), df, weights = w,
           control = el_control(tol = 1e-08, th = 1e+10, nthreads = 1)))
  wlfit2 <- suppressWarnings(
    el_glm(y ~ ., family = binomial(link = "logit"), df, weights = w,
           control = el_control(tol = 1e-08, th = 1e+10)))
  expect_equal(wlfit@optim, wlfit2@optim)
  expect_equal(wlfit@parTests, wlfit2@parTests)

  pfit <- el_glm(y ~ ., family = binomial(link = "probit"), df,
                control = el_control(tol = 1e-08, th = 1e+10, nthreads = 1))
  pfit2 <- el_glm(y ~ ., family = binomial(link = "probit"), df,
                 control = el_control(tol = 1e-08, th = 1e+10))
  expect_equal(pfit@optim, pfit2@optim)
  expect_equal(pfit@parTests, pfit2@parTests)
  wpfit <- suppressWarnings(
    el_glm(y ~ ., family = binomial(link = "probit"), df, weights = w,
           control = el_control(tol = 1e-08, th = 1e+10, nthreads = 1)))
  wpfit2 <- suppressWarnings(
    el_glm(y ~ ., family = binomial(link = "probit"), df, weights = w,
           control = el_control(tol = 1e-08, th = 1e+10)))
  expect_equal(wpfit@optim, wpfit2@optim)
  expect_equal(wpfit@parTests, wpfit2@parTests)
})
