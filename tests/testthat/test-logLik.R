test_that("logLik at maximum EL estimates", {
  data("airquality")
  n <- nrow(airquality)
  fit <- el_lm(Wind ~ Temp, data = airquality)
  loglik <- suppressWarnings(logLik(fit, REML = T))
  expect_output(show(loglik))
  expect_output(print(loglik))
  expect_equal(loglik@logLik, -n * log(n))
})

test_that("empty model", {
  data("mtcars")
  fit <- el_lm(mpg ~ 0, mtcars)
  expect_error(logLik(fit))
})
