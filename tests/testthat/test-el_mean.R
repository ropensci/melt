test_that("invalid 'x", {
  expect_error(el_mean(c(1, Inf), par = 0))
  expect_error(el_mean(10, 0))
  expect_error(el_mean(matrix(c(1, 1, 2, 2), ncol = 2), par = c(0, 0)))
  expect_error(el_mean(matrix(c(1, 1, 2, 2), ncol = 2), par = c(0, 0),
    weights = c(1, 2)
  ))
})

test_that("invalid 'par", {
  x <- sleep$extra
  expect_error(el_mean(x, par = NA))
  expect_error(el_mean(x, par = NULL))
  expect_error(el_mean(x, par = Inf))
  expect_error(el_mean(x, par = NaN))
  expect_error(el_mean(x, par = c(0, 0)))
})

test_that("convergence check", {
  x <- women$weight
  grid <- seq(120, 160, length.out = 1000)
  conv <- function(par) {
    el_mean(x, par)@optim$convergence
  }
  expect_true(all(vapply(grid, conv, FUN.VALUE = logical(1))))
})

test_that("probabilities add up to 1", {
  x <- women$height
  w <- women$weight
  fit <- el_mean(x, par = 65)
  expect_output(print(fit))
  expect_equal(sum(exp(fit@logp)), 1, tolerance = 1e-07)
  fit2 <- el_mean(x, par = 65, weights = w)
  expect_equal(sum(exp(fit2@logp)), 1, tolerance = 1e-7)
})

test_that("identical weights == no weights", {
  x <- women$height
  fit <- el_mean(x, par = 65)
  fit2 <- el_mean(x, par = 65, weights = rep(0.5, length(x)))
  fit2@weights <- numeric()
  expect_equal(fit, fit2)
})

test_that("conversion between loglik and loglr", {
  x <- women$height
  n <- length(x)
  fit <- el_mean(x, par = 65)
  expect_equal(fit@logl + n * log(n), fit@loglr)
  fit2 <- el_mean(x, par = 65, weights = women$weight)
  w <- weights(fit2)
  expect_equal(fit2@logl + sum(w * (log(n) - log(w))), fit2@loglr)
})

test_that("identical weights == no weights", {
  x <- women$height
  w <- rep(mean(women$weight), length(x))
  par <- 60
  fit <- el_mean(x, par)
  fit2 <- el_mean(x, par, weights = w)
  expect_equal(fit@optim, fit2@optim)
})

# test_that("same results from repeated executions", {
#   x <- women$height
#   w <- women$weight
#   fit <- el_mean(x, par = 65)
#   fit2 <- el_mean(x, par = 65)
#   expect_equal(fit, fit2)
#
#   wfit <- el_mean(x, par = 65, weights = w)
#   wfit2 <- el_mean(x, par = 65, weights = w)
#   expect_equal(wfit, wfit2)
#
#   lhs <- 4
#   elt <- elt(fit, lhs = lhs)
#   elt2 <- elt(fit2, lhs = lhs)
#   expect_equal(elt, elt2)
#
#   elt3 <- elt(wfit, lhs = lhs)
#   elt4 <- elt(wfit2, lhs = lhs)
#   expect_equal(elt3, elt4)
# })
