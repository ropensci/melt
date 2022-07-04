test_that("Invalid `formula`", {
  fit <- el_lm(mpg ~ 0, data = mtcars)
  expect_error(logLik(fit))
})

test_that("`logLik()` at the maximum empirical likelihood estimate.", {
  n <- nrow(airquality)
  fit <- el_lm(Wind ~ Temp, data = airquality)
  loglik <- suppressWarnings(logLik(fit, REML = TRUE))
  expect_output(show(loglik))
  expect_output(print(loglik))
  expect_equal(loglik@logLik, -n * log(n))
})
