test_that("Invalid `x`.", {
  expect_error(el_mean(c(1, Inf), par = 0))
  expect_error(el_mean(10, 0))
  expect_error(el_mean(c(1, 1), par = 1))
  expect_error(el_mean(matrix(c(1, 1, 2, 2), ncol = 2), par = c(0, 0)))
  expect_error(el_mean(matrix(c(1, 1, 2, 2), ncol = 2),
    par = c(0, 0),
    weights = c(1, 2)
  ))
})

test_that("Invalid `par`.", {
  expect_error(el_mean(precip, par = NA))
  expect_error(el_mean(precip, par = NULL))
  expect_error(el_mean(precip, par = Inf))
  expect_error(el_mean(precip, par = NaN))
  expect_error(el_mean(precip, par = c(0, 0)))
})

test_that("Convergence check.", {
  x <- women$weight
  grid <- seq(10, 60, length.out = 100)
  expect_true(all(vapply(grid, function(par) {
    conv(el_mean(precip, par))
  },
  FUN.VALUE = logical(1)
  )))
})

test_that("Probabilities add up to 1.", {
  x <- women$height
  w <- women$weight
  fit <- el_mean(x, par = 60)
  expect_equal(sum(exp(logProb(fit))), 1, tolerance = 1e-07)
  fit2 <- el_mean(x, par = 60, weights = w)
  expect_equal(sum(exp(logProb(fit2))), 1, tolerance = 1e-07)
})

test_that("Identical weights means no weights.", {
  fit <- el_mean(precip, par = 60)
  fit2 <- el_mean(precip, par = 60, weights = rep(0.5, length(precip)))
  fit2@weights <- numeric()
  expect_equal(getOptim(fit), getOptim(fit2))
})

test_that("Conversion between `logl` and `loglr`.", {
  x <- women$height
  n <- length(x)
  fit <- el_mean(x, par = 60)
  expect_equal(logL(fit) + n * log(n), logLR(fit))
  fit2 <- el_mean(x, par = 60, weights = women$weight)
  w <- weights(fit2)
  expect_equal(logL(fit2) + sum(w * (log(n) - log(w))), logLR(fit2))
})

test_that("`verbose` == TRUE in `el_control()`.", {
  expect_message(el_mean(precip,
    par = 60, control = el_control(verbose = TRUE)
  ))
})

test_that("`conv()` methods.", {
  fit <- el_mean(precip, par = 60)
  fit2 <- el_mean(precip, par = 0)
  expect_true(conv(fit))
  expect_false(conv(fit2))
})

test_that("Larger `tol_l` decreases iterations for convergence.", {
  fit <- el_mean(precip, par = 60, control = el_control(tol_l = 1e-08))
  fit2 <- el_mean(precip, par = 60, control = el_control(tol_l = 1e-02))
  expect_gte(getOptim(fit)$iterations, getOptim(fit2)$iterations)
})

test_that("Noise susceptibility tests.", {
  fit <- el_mean(precip, par = 60)
  fit2 <- el_mean(precip, par = 60 + .Machine$double.eps)
  expect_equal(getOptim(fit), getOptim(fit2))
})

test_that("Convex hull constraint violated.", {
  x <- women$weight
  grid <- seq(70, 100, length.out = 100)
  conv <- function(par) {
    getOptim(el_mean(precip, par))$convergence
  }
  expect_false(any(vapply(grid, conv, FUN.VALUE = logical(1))))
})

test_that("`print()` method.", {
  fit <- el_mean(precip, par = 60)
  expect_output(show(fit))
  expect_output(print(fit))
  out <- summary(fit)
  expect_output(print(out))
  expect_output(show(out))
  fit@statistic <- numeric()
  expect_output(print(fit))
  fit2 <- el_mean(women$height, par = 60, weights = women$weight)
  out2 <- summary(fit2)
  expect_output(print(out2))
  expect_output(show(out2))
})
