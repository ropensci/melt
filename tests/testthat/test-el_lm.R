test_that("invalid 'formula'", {
  expect_error(el_lm(cbind(mtcars$mpg, mtcars$cyl) ~ 1, data = mtcars))
})

test_that("invalid 'data'", {
  x <- mtcars$disp
  x2 <- x
  y <- mtcars$mpg
  df <- data.frame(y, x, x2)
  w <- mtcars$carb
  expect_error(el_lm(y ~ x + x2, data = df))
  expect_error(el_lm(y ~ x + x2, data = df, weights = w))
})

test_that("invalid 'weights'", {
  w <- women$weight
  w[1] <- -1
  expect_error(el_lm(height ~ weight, data = women, weights = w))
  w[1] <- "error"
  expect_error(el_lm(height ~ weight, data = women, weights = w))
})

test_that("invalid 'control'", {
  expect_error(el_lm(mpg ~ ., data = mtcars, control = list()))
})

test_that("probabilities add up to 1", {
  x <- women$weight
  y <- 10 + 0.001 * x + rep(c(1, 0.2, 0.5, 2, -1.2), times = 3)
  df <- data.frame(x, y)
  fit <- el_lm(y ~ x, data = df)
  expect_equal(sum(exp(fit@logp)), 1)
  fit2 <- el_lm(y ~ x, data = df, weights = women$height)
  expect_equal(sum(exp(fit2@logp)), 1)
})

test_that("conversion between loglik and loglr", {
  fit <- el_lm(eruptions ~ waiting, data = faithful)
  n <- nrow(faithful)
  expect_equal(fit@logl + n * log(n), fit@loglr)
  fit2 <- el_lm(eruptions ~ waiting,
    data = faithful,
    weights = faithful$waiting
  )
  w <- weights(fit2)
  expect_equal(fit2@logl + sum(w * (log(n) - log(w))), fit2@loglr)
})

test_that("empty model", {
  fit <- el_lm(eruptions ~ 0, faithful)
  expect_output(print(summary(fit)))
})

test_that("no intercept", {
  fit <- el_lm(mpg ~ -1 + ., mtcars)
  expect_s4_class(fit, "LM")
})

test_that("print method", {
  fit <- el_lm(mpg ~ disp + hp, data = mtcars)
  expect_output(show(fit))
  expect_output(print(fit))
  out <- summary(fit)
  expect_output(print(out))
  expect_output(show(out))
  out@aliased <- c(TRUE, TRUE, TRUE)
  expect_output(print(out))
  fit@statistic <- numeric()
  expect_output(print(fit))
  df2 <- mtcars
  df2[1, 1] <- NA
  fit2 <- el_lm(mpg ~ disp + hp, data = df2)
  expect_output(print(summary(fit2)))
})







# test_that("print method", {
#   n <- 100
#   x <- rnorm(n)
#   x2 <- rnorm(n)
#   y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
#   df <- data.frame(y, x, x2)
#   optcfg <- el_control(tol = 1e-08, th = 1e+10)
#   fit <- el_lm(y ~ x + x2, df, control = optcfg)
#   expect_output(show(fit))
#   expect_output(print(fit))
#   out <- summary(fit)
#   expect_output(print(out))
#   expect_output(show(out))
#   fit@statistic <- numeric()
#   expect_output(print(fit))
#
#   df2 <- df
#   df2[1, 1] <- NA
#   fit2 <- el_lm(y ~ x + x2, df2, control = optcfg)
#   expect_output(print(summary(fit2)))
#
#   out@aliased <- c(TRUE, TRUE, TRUE)
#   expect_output(print(out))
# })

# test_that("same results with parallel computing", {
#   skip_on_os("windows", arch = "i386")
#   n <- 400
#   p <- 15
#   b <- rnorm(p)
#   x <- matrix(rnorm(n * p), ncol = p)
#   y <- 1 + x %*% as.vector(b) + rnorm(n)
#   df <- data.frame(y, x)
#   fit <- el_lm(y ~ ., df, control = el_control(th = 1e+10, nthreads = 1))
#   fit2 <- el_lm(y ~ ., df, control = el_control(th = 1e+10))
#   expect_equal(fit@optim, fit2@optim)
#   expect_equal(fit@parTests, fit2@parTests)
#
#   w <- 1 + runif(n, min = -0.5, max = 0.5)
#   wfit <- el_lm(y ~ ., df,
#     weights = w,
#     control = el_control(th = 1e+10, nthreads = 1)
#   )
#   wfit2 <- el_lm(y ~ ., df,
#     weights = w,
#     control = el_control(th = 1e+10)
#   )
#   expect_equal(wfit@optim, wfit2@optim)
#   expect_equal(wfit@parTests, wfit2@parTests)
# })
