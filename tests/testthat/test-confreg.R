test_that("normal", {
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.4 * x - 0.2 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  cr <- confreg(fit, parm = c(1, 2))
  pdf(NULL)
  plot(cr)
  expect_s4_class(cr, "ConfregEL")
})

test_that("empty model", {
  data("mtcars")
  fit <- el_lm(mpg ~ 0, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2)))
})

test_that("one parameter", {
  data("women")
  fit <- el_mean(women$height, 67)
  expect_error(confreg(fit, parm = c(1, 2)))
})

test_that("invalid 'parm'", {
  data("mtcars")
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2, 3)))
  expect_error(confreg(fit, parm = c(1, 1)))
  expect_error(confreg(fit, parm = c("error", "error2")))
  expect_error(confreg(fit, parm = c(NaN, NA)))
})

test_that("invalid 'level'", {
  data("mtcars")
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_s4_class(confreg(fit, parm = c(1, 2), level = 0), "ConfregEL")
  expect_error(confreg(fit, parm = c(1, 2), level = 1))
})

test_that("invalid 'cv'", {
  data("mtcars")
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2), cv = 1e+10))
})

test_that("invalid 'npoints'", {
  data("mtcars")
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2), npoints = -10))
})
