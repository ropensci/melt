test_that("invalid 'g'", {
  expect_error(el_eval(matrix(c(1, 1), ncol = 2)))
  expect_error(el_eval(matrix(c(1, 1, 2, NA), ncol = 2)))
  expect_error(el_eval(matrix(c(1, 1, 2, 2), ncol = 2)))
  expect_error(el_eval(matrix(c(1, 1, 2, 2), ncol = 2), weights = c(1, 2)))
})

test_that("invalid 'control'", {
  expect_error(el_eval(women$height - 67, control = list(maxit = 200L)))
})

test_that("convergence check", {
  x <- women$weight
  grid <- seq(120, 160, length.out = 1000)
  conv <- function(par) {
    el_eval(x - par)$optim$convergence
  }
  expect_true(all(vapply(grid, conv, FUN.VALUE = logical(1))))
})

test_that("probabilities add up to 1", {
  x <- women$height
  w <- women$weight
  fit <- el_eval(x - 60)
  expect_equal(sum(exp(fit$logp)), 1, tolerance = 1e-07)
  fit2 <- el_eval(x - 60, weights = w)
  expect_equal(sum(exp(fit2$logp)), 1, tolerance = 1e-07)
})

test_that("conversion between loglik and loglr", {
  x <- women$height
  n <- length(x)
  w <- women$weight
  fit <- el_eval(x - 60)
  expect_equal(fit$logl + n * log(n), fit$loglr, tolerance = 1e-07)
  fit2 <- el_eval(x - 60, weights = w)
  w <- weights(fit2)
  expect_equal(fit2$logl + sum(w * (log(n) - log(w))), fit2$loglr,
    tolerance = 1e-07
  )
})
