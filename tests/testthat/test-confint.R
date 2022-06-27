test_that("invalid 'parm'", {
  data("women")
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(women$height, 67, control = optcfg)
  expect_error(confint(fit, parm = NA))
  expect_error(confint(fit, parm = NULL))
  expect_error(confint(fit, parm = NaN))
  expect_error(confint(fit, parm = Inf))
})

test_that("invalid 'level'", {
  data("women")
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(women$height, 67, control = optcfg)
  expect_error(confint(fit, level = -1))
  expect_error(confint(fit, level = 2))
  expect_error(confint(fit, level = Inf))
  expect_error(confint(fit, level = c(0, 0)))
})

test_that("invalid 'control'", {
  data("women")
  fit <- el_mean(women$height, 67)
  expect_error(confint(fit, control = list(maxit = 67)))
})

test_that("'level' == 1", {
  data("women")
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(women$height, 67, control = optcfg)
  ci <- confint(fit, level = 1)
  expect_equal(ci[1], -Inf)
  expect_equal(ci[2], Inf)
})

test_that("'level' == 0", {
  data("women")
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(women$height, 67, control = optcfg)
  ci <- confint(fit, level = 0)
  expect_equal(ci[2] - ci[1], 0)
})

test_that("empty model", {
  data("mtcars")
  fit <- el_lm(mpg ~ -1, mtcars)
  ci <- confint(fit)
  expect_equal(nrow(ci), 0)
})

test_that("'nthreads' == 1", {
  data("mtcars")
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

test_that("unnamed coefficients", {
  data("faithful")
  fit <- el_mean(as.matrix(faithful), par = c(4, 60))
  names(fit@coefficients) <- NULL
  expect_type(confint(fit, parm = c(1, 2)), "double")
  expect_type(confint(fit), "double")
})

test_that("character specification for 'parm'", {
  data("faithful")
  fit <- el_mean(as.matrix(faithful), par = c(4, 60))
  expect_type(confint(fit, parm = c("eruptions2", "eruptions3")), "double")
  expect_type(confint(fit, parm = c("eruptions", "eruptions2")), "double")
})
