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
  x <- sleep$extra
  expect_error(el_mean(x, par = NA))
  expect_error(el_mean(x, par = NULL))
  expect_error(el_mean(x, par = Inf))
  expect_error(el_mean(x, par = NaN))
  expect_error(el_mean(x, par = c(0, 0)))
})

test_that("Convergence check.", {
  x <- women$weight
  grid <- seq(120, 160, length.out = 1000)
  expect_true(all(vapply(grid, function(par) {
    conv(el_mean(x, par))
  },
  FUN.VALUE = logical(1)
  )))
})

test_that("Probabilities add up to 1.", {
  x <- women$height
  w <- women$weight
  fit <- el_mean(x, par = 60)
  expect_output(print(fit))
  expect_equal(sum(exp(fit@logp)), 1, tolerance = 1e-07)
  fit2 <- el_mean(x, par = 60, weights = w)
  expect_equal(sum(exp(fit2@logp)), 1, tolerance = 1e-07)
})

test_that("Identical weights means no weights.", {
  x <- women$height
  fit <- el_mean(x, par = 60)
  fit2 <- el_mean(x, par = 60, weights = rep(0.5, length(x)))
  fit2@weights <- numeric()
  expect_equal(fit, fit2)
})

test_that("Conversion between `loglik` and `loglr`.", {
  x <- women$height
  n <- length(x)
  fit <- el_mean(x, par = 60)
  expect_equal(fit@logl + n * log(n), logLR(fit))
  fit2 <- el_mean(x, par = 60, weights = women$weight)
  w <- weights(fit2)
  expect_equal(fit2@logl + sum(w * (log(n) - log(w))), logLR(fit2))
})

test_that("`verbose` == TRUE in `el_control()`.", {
  x <- women$height
  expect_message(el_mean(x, par = 60, control = el_control(verbose = TRUE)))
})

test_that("`conv()` methods.", {
  x <- women$height
  fit <- el_mean(x, par = 60)
  expect_true(conv(fit))
})

#' @srrstats {G5.7} Larger `tol_l` decreases the number of iterations for
#'   convergence in `el_mean()`.
test_that("Larger `tol_l` decreases iterations for convergence.", {
  x <- women$height
  fit <- el_mean(x, par = 60, control = el_control(tol_l = 1e-08))
  fit2 <- el_mean(x, par = 60, control = el_control(tol_l = 1e-02))
  expect_gte(fit@optim$iterations, fit2@optim$iterations)
})

#' @srrstats {G5.9, G5.9a} Adding trivial noise does not change the overall
#'   optimization results.
test_that("Noise susceptibility tests.", {
  x <- women$height
  fit <- el_mean(x, par = 60)
  fit2 <- el_mean(x, par = 60 + .Machine$double.eps)
  expect_equal(fit@optim, fit2@optim)
})

#' @srrstats {RE1.4} Violation of the convex hull constraint is tested.
test_that("Convex hull constraint violated.", {
  x <- women$weight
  grid <- seq(10, 50, length.out = 1000)
  conv <- function(par) {
    el_mean(x, par)@optim$convergence
  }
  expect_false(any(vapply(grid, conv, FUN.VALUE = logical(1))))
})
