test_that("invalid 'level'", {
  data("women")
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(women$height, 67, control = optcfg)
  expect_error(confint(fit, level = -1))
  expect_error(confint(fit, level = 2))
  expect_error(confint(fit, level = Inf))
  expect_error(confint(fit, level = c(0, 0)))
})

test_that("invalid 'parm'", {
  data("women")
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(women$height, 67, control = optcfg)
  expect_error(confint(fit, parm = NA))
  expect_error(confint(fit, parm = NULL))
  expect_error(confint(fit, parm = NaN))
  expect_error(confint(fit, parm = Inf))
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

test_that("no effect of nthreads", {
  data("mtcars")
  fit <- el_lm(mpg ~ disp + hp + wt + qsec, mtcars)
  parm <- rep(c(2, 5, 1), times = 10)
  ci1 <- confint(fit, parm = parm, control = el_control(nthreads = 1L))
  ci2 <- confint(fit, parm = parm)
  expect_equal(ci1, ci2)
})
