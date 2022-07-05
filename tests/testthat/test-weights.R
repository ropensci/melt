test_that("Invalid `weights`.", {
  x <- airquality$Wind
  w <- airquality$Day
  expect_error(el_mean(x, par = 10, weights = w[-1]))
  w[1] <- -10
  expect_error(el_mean(x, par = 10, weights = w))
  w[1] <- Inf
  expect_error(el_mean(x, par = 10, weights = w))
  expect_error(el_mean(x, par = 10, weights = rep("w", length(x))))
})

test_that("Re-scaled weights.", {
  x <- airquality$Wind
  w <- airquality$Day
  fit <- el_mean(x, par = 10, control = el_control(step = 1, maxit_l = 100))
  expect_warning(weights(fit, "extra arguments"))
  expect_identical(weights(fit), NULL)
  wfit <- el_mean(x, par = 10, weights = w)
  expect_equal(sum(weights(wfit)), length(x))
  expect_output(print(wfit))
})
