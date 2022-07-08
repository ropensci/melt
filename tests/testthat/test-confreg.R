test_that("invalid 'object'", {
  fit <- el_lm(mpg ~ 0, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2)))
})

test_that("invalid 'parm'", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2, 3)))
  expect_error(confreg(fit, parm = c(1, 1)))
  expect_error(confreg(fit, parm = c("error", "error2")))
  expect_error(confreg(fit, parm = c(NaN, NA)))
  names(fit@coefficients) <- NULL
  expect_error(confreg(fit, parm = c("error", "error2")))
  fit2 <- el_mean(women$height, 67)
  expect_error(confreg(fit2, parm = c(1, 2)))
})

test_that("invalid 'level'", {
  fit <- el_lm(mpg ~ hp + disp, data = mtcars)
  expect_s4_class(confreg(fit, parm = c(1, 2), level = 0), "ConfregEL")
  expect_error(confreg(fit, parm = c(1, 2), level = 1))
})

test_that("invalid 'cv'", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2), cv = 1e+10))
  expect_error(confreg(fit,
                       parm = c(1, 2), cv = 1e+10,
                       control = el_control(th = 100)
  ))
})

test_that("invalid 'npoints'", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2), npoints = -10))
})

test_that("invalid 'control'", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2), control = list(maxit = 10)))
})

test_that("plot method", {
  fit <- el_lm(mpg ~ disp + hp + wt, data = mtcars)
  cr <- confreg(fit)
  cr2 <- confreg(fit, parm = c(1, 2), cv = 6)
  pdf(NULL)
  expect_invisible(plot(cr2))
})


