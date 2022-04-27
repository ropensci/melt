test_that("invalid 'weights'", {
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  w <- rep(runif(1), length(x) - 1L)
  expect_error(el_mean(par, x, weights = w))
  w <- rep(runif(1), length(x))
  w[1] <- -10
  expect_error(el_mean(par, x, weights = w))
  w[1] <- Inf
  expect_error(el_mean(par, x, weights = w))
  expect_error(el_mean(par, x, weights = rep("w", length(x))))
})
