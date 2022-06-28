test_that("invalid 'object'", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  fit@data <- matrix(NA_real_, nrow = 0L, ncol = 0L)
  expect_error(eld(fit))
})

test_that("invalid 'control'", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(eld(fit, control = list(maxit = 20L)))
})

test_that("weights", {
  fit <- el_lm(extra ~ ID, data = sleep, weights = as.numeric(group))
  expect_s4_class(eld(fit), "ELD")
})

test_that("mean", {
  fit <- el_mean(women$height, par = 67)
  fit2 <- el_mean(women$height, par = 67, weights = women$weight)
  expect_s4_class(eld(fit), "ELD")
  expect_s4_class(eld(fit2), "ELD")
})

test_that("plot method", {
  fit <- el_lm(mpg ~ disp + hp + wt, data = mtcars)
  eld <- eld(fit)
  pdf(NULL)
  expect_invisible(plot(eld))
})
