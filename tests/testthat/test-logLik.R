test_that("Invalid `formula`.", {
  fit <- el_lm(mpg ~ 0, data = mtcars)
  expect_error(logLik(fit))
})

test_that("`logLik()` at the maximum empirical likelihood estimate.", {
  n <- nrow(airquality)
  fit <- el_lm(Wind ~ Temp, data = airquality)
  logLik <- suppressWarnings(logLik(fit, REML = TRUE))
  expect_output(show(logLik))
  expect_output(print(logLik))
  expect_equal(getDataPart(logLik), -n * log(n))
})
