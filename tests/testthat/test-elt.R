test_that("Invalid `object`.", {
  fit <- el_eval(mtcars$mpg)
  expect_error(elt(fit, lhs = 1))
  fit2 <- el_lm(mpg ~ 0, data = mtcars)
  expect_error(elt(fit2, lhs = 1))
})

test_that("Invalid `rhs` and `lhs`.", {
  fit <- el_lm(mpg ~ disp + hp, data = mtcars)
  fit2 <- el_lm(mpg ~ disp + hp, data = mtcars, weights = mtcars$wt)
  lhs <- matrix(c(0, 1, -1), nrow = 1)
  expect_error(elt(fit))
  expect_error(elt(fit, rhs = c(NA, 0)))
  expect_error(elt(fit, rhs = c(1, 0), lhs = lhs))
  expect_error(elt(fit, rhs = matrix(c(1, 0, 0), ncol = 3), lhs = lhs))
  expect_error(elt(fit, rhs = matrix(c(1, 0, 0), ncol = 1), lhs = lhs))
  expect_error(elt(fit, rhs = matrix(c("error"), ncol = 1), lhs = lhs))
  expect_error(elt(fit, lhs = matrix(c(1, 0, 0, 0, NA, 1), nrow = 2)))
  expect_error(elt(fit2, lhs = matrix(c(1, 0, 0, 0, NA, 1), nrow = 2)))
  expect_error(elt(fit, lhs = matrix(rnorm(4), ncol = 2)))
  expect_error(elt(fit2, lhs = matrix(rnorm(4), ncol = 2)))
  expect_error(elt(fit, lhs = matrix(rnorm(12), ncol = 3)))
  expect_error(elt(fit2, lhs = matrix(rnorm(12), ncol = 3)))
  expect_error(elt(fit, lhs = matrix(c(1, 1, 0, 0, 0, 0), nrow = 2)))
  expect_error(elt(fit2, lhs = matrix(c(1, 1, 0, 0, 0, 0), nrow = 2)))
})

test_that("Invalid `calibrate`.", {
  fit <- el_lm(mpg ~ wt, data = mtcars)
  expect_error(elt(fit, rhs = c(1, 1), calibrate = c(1, 2)))
  expect_error(elt(fit, rhs = c(1, 1), calibrate = "error"))
  expect_error(elt(fit, rhs = c(1, 1), calibrate = "f"))
})

test_that("Invalid `control`.", {
  fit <- el_mean(sleep$extra, par = 0)
  expect_error(elt(fit, lhs = 1, control = list(maxit = 200L)))
})

test_that("When elt == evaluation.", {
  x <- sleep$extra
  fit <- el_mean(x, par = 1.2)
  fit2 <- elt(fit, rhs = 1.2)
  expect_equal(fit@optim$lambda, fit2@optim$lambda)
  wfit <- el_mean(x, par = 1.2, weights = as.numeric(sleep$group))
  wfit2 <- elt(wfit, rhs = 1.2)
  expect_equal(wfit@optim$lambda, wfit2@optim$lambda)
  fit3 <- el_lm(mpg ~ disp + hp, data = mtcars)
  lhs <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 2)
  rhs <- c(0, 0)
  fit4 <- elt(fit3, rhs = rhs, lhs = lhs)
  expect_output(print(fit3))
  expect_equal(fit3@optim$lambda, fit4@optim$lambda)
  wfit3 <- el_lm(mpg ~ disp + hp, data = mtcars, weights = wt)
  lhs <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 2)
  rhs <- c(0, 0)
  wfit4 <- elt(wfit3, rhs = rhs, lhs = lhs)
  expect_equal(wfit3@optim$lambda, wfit4@optim$lambda)
})

test_that("`conv()` method and calibration.", {
  fit <- el_mean(women$height, par = 65)
  expect_true(conv(elt(fit, rhs = 67, calibrate = "f")))
  expect_s4_class(elt(fit, rhs = 67, calibrate = "boot"), "ELT")
})

test_that("Vector `lhs`.", {
  fit <- el_lm(mpg ~ disp + hp, data = mtcars)
  expect_error(elt(fit, lhs = c(1, -1, 0, 0)))
  out <- elt(fit, lhs = c(1, -1, 0))
  expect_s4_class(out, "ELT")
  expect_output(show(out))
  expect_output(print(out))
  fit2 <- el_mean(faithful, par = c(4, 70), weights = faithful$waiting)
  out2 <- elt(fit2, lhs = c(17, -1), rhs = -3)
  expect_s4_class(out, "ELT")
  expect_output(show(out))
  expect_output(print(out))
})

test_that("Matrix `rhs`.", {
  fit <- el_lm(mpg ~ disp + hp, data = mtcars)
  rhs <- matrix(c(0, 1, 2), ncol = 1)
  out <- suppressMessages(elt(fit, rhs = rhs))
  expect_s4_class(out, "ELT")
})

test_that("Missing `object`.", {
  expect_null(elt(rhs = 1))
})

test_that("`SD` class.", {
  x <- women$height
  fit <- el_sd(x, mean = 65, sd = 4)
  fit@npar <- 0L
  expect_error(elt(fit, rhs = 1))
  fit@npar <- 1L
  fit@data <- NULL
  expect_error(elt(fit, rhs = 1))
  fit <- el_sd(x, mean = 65, sd = 4)
  expect_error(elt(fit, rhs = 1, control = list()))
  expect_error(elt(fit, rhs = 1, calibrate = "f"))
  expect_error(elt(fit, rhs = "error"))
  expect_error(elt(fit, rhs = c(1, 2)))
  expect_error(elt(fit, rhs = -1))
  expect_error(elt(fit, rhs = rhs, lhs = lhs, calibrate = "boot"))
  expect_error(elt(fit, rhs = rhs, lhs = lhs, calibrate = "f"))
  expect_error(elt(fit, rhs = -1, lhs = 1, calibrate = "f"))
  out <- elt(fit, rhs = 1)
  out2 <- elt(fit, rhs = 1, lhs = 2)
  expect_s4_class(out, "ELT")
  expect_s4_class(out2, "ELT")
})
