test_that("invalid 'weights'", {
  data("airquality")
  x <- airquality$Wind
  w <- airquality$Day
  par <- 10
  expect_error(el_mean(x, par, weights = w[-1]))
  w[1] <- -10
  expect_error(el_mean(x, par, weights = w))
  w[1] <- Inf
  expect_error(el_mean(x, par, weights = w))
  expect_error(el_mean(x, par, weights = rep("w", length(x))))
})

test_that("weights() function", {
  data("airquality")
  x <- airquality$Wind
  w <- airquality$Day
  n <- length(x)
  fit <- el_mean(x, par = 10, control = el_control(step = 1, maxit_l = 100))
  expect_warning(weights(fit, "extra arguments"))
  expect_identical(weights(fit), NULL)
  fit2 <- el_mean(x, par = 10, weights = w)
  expect_equal(sum(weights(fit2)), length(x))
})
