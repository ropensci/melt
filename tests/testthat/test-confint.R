test_that("Invalid `parm`.", {
  fit <- el_mean(women$height, par = 67)
  fit2 <- el_sd(women$height, mean = 65, sd = 5)
  expect_error(confint(fit, parm = NA))
  expect_error(confint(fit, parm = NULL))
  expect_error(confint(fit, parm = NaN))
  expect_error(confint(fit, parm = Inf))
  expect_error(confint(fit2, parm = NA))
  expect_error(confint(fit2, parm = NULL))
  expect_error(confint(fit2, parm = NaN))
  expect_error(confint(fit2, parm = Inf))
})

test_that("Invalid `level`.", {
  fit <- el_mean(women$height, par = 67)
  fit2 <- el_sd(women$height, mean = 65, sd = 5)
  expect_error(confint(fit, level = -1))
  expect_error(confint(fit, level = 2))
  expect_error(confint(fit, level = Inf))
  expect_error(confint(fit, level = c(0, 0)))
  expect_error(confint(fit2, level = -1))
  expect_error(confint(fit2, level = 2))
  expect_error(confint(fit2, level = Inf))
  expect_error(confint(fit2, level = c(0, 0)))
})

test_that("Invalid `control`.", {
  fit <- el_mean(women$height, par = 67)
  fit2 <- el_sd(women$height, mean = 65, sd = 5)
  expect_error(confint(fit, control = list(maxit = 67)))
  expect_error(confint(fit2, control = list(maxit = 67)))
})

test_that("`level` == 1.", {
  fit <- el_mean(women$height, par = 67)
  fit2 <- el_sd(women$height, mean = 65, sd = 5)
  ci <- confint(fit, level = 1)
  ci2 <- confint(fit2, level = 1)
  expect_equal(ci[1], -Inf)
  expect_equal(ci[2], Inf)
  expect_equal(ci2[1], 0)
  expect_equal(ci2[2], Inf)
})

test_that("`level` == 0.", {
  fit <- el_mean(women$height, par = 67)
  fit2 <- el_sd(women$height, mean = 65, sd = 5)
  ci <- confint(fit, level = 0)
  ci2 <- confint(fit2, level = 0)
  expect_equal(ci[2] - ci[1], 0)
  expect_equal(ci2[2] - ci2[1], 0)
})

test_that("Empty model.", {
  fit <- el_lm(mpg ~ -1, data = mtcars)
  ci <- confint(fit)
  expect_equal(nrow(ci), 0)
})

test_that("`nthreads` == 1.", {
  mpg <- mtcars$mpg
  disp <- mtcars$disp
  hp <- mtcars$hp
  wt <- mtcars$wt
  qsec <- mtcars$qsec
  fit <- el_lm(mpg ~ disp + hp + wt + qsec)
  parm <- rep(c(2, 5, 1), times = 10)
  ci1 <- confint(fit, parm = parm, control = el_control(nthreads = 1L))
  ci2 <- confint(fit, parm = parm)
  expect_equal(ci1, ci2)
})

test_that("Unnamed coefficients.", {
  fit <- el_mean(as.matrix(faithful), par = c(4, 60))
  fit2 <- el_sd(women$height, mean = 65, sd = 5)
  names(fit@coefficients) <- NULL
  expect_type(confint(fit, parm = c(1, 2)), "double")
  expect_type(confint(fit), "double")
  expect_type(confint(fit2, parm = c(1, 2)), "double")
  expect_type(confint(fit2), "double")
})

test_that("Character specification for `parm`.", {
  fit <- el_mean(faithful, par = c(4, 60))
  fit2 <- el_sd(women$height, mean = 35, sd = 13)
  fit3 <- el_sd(women$height, mean = 35, sd = 13, weights = women$weight)
  expect_type(confint(fit, parm = c("eruptions2", "eruptions3")), "double")
  expect_type(confint(fit, parm = c("eruptions", "eruptions2")), "double")
  expect_type(confint(fit2, parm = "1"), "double")
  expect_type(confint(fit2, parm = c("e", "f")), "double")
  expect_type(confint(fit2, parm = c("1", "f")), "double")
  expect_type(confint(fit3, parm = "1"), "double")
  names(fit2@coefficients) <- "it"
  expect_type(confint(fit2, parm = 1), "double")
})
