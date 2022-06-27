test_that("invalid 'object", {
  n <- nrow(mtcars)
  fit <- el_lm(mpg ~ hp, mtcars)
  expect_error(eld(fit, control = list(maxit = 20L)))
  eld <- eld(fit)
  pdf(NULL)
  plot(eld)
  expect_length(eld@eld, n)
  lhs <- matrix(c(1, -1), nrow = 1)
  fit2 <- elt(fit, lhs = lhs)
  expect_error(eld(mtcars))
  fit@data <- matrix(NA_real_, nrow = 0L, ncol = 0L)
  expect_error(eld(fit))
  expect_error(eld(fit2))
})

test_that("weights", {
  fit <- el_lm(extra ~ ID, sleep, weights = as.numeric(group))
  expect_s4_class(eld(fit), "ELD")
})

test_that("mean", {
  fit <- el_mean(women$height, par = 67)
  fit2 <- el_mean(women$height, par = 67, weights = women$weight)
  expect_s4_class(eld(fit), "ELD")
  expect_s4_class(eld(fit2), "ELD")
})
