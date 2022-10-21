test_that("Invalid `parm`.", {
  fit <- el_mean(precip, par = 60)
  expect_error(confint(fit, parm = NA))
  expect_error(confint(fit, parm = NULL))
  expect_error(confint(fit, parm = NaN))
  expect_error(confint(fit, parm = Inf))
})

test_that("Invalid `level`.", {
  fit <- el_mean(precip, par = 60)
  expect_error(confint(fit, level = -1))
  expect_error(confint(fit, level = 2))
  expect_error(confint(fit, level = Inf))
  expect_error(confint(fit, level = c(0, 0)))
})

test_that("Invalid `control`.", {
  fit <- el_mean(precip, par = 60)
  expect_error(confint(fit, control = list(maxit = 67)))
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
  expect_equal(nrow(ci), 0L)
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
  fit <- el_sd(women$height, mean = 65, sd = 5)
  expect_type(confint(fit, parm = c(1, 2)), "double")
  expect_type(confint(fit), "double")
})

test_that("Character specification for `parm`.", {
  fit <- el_sd(women$height, mean = 35, sd = 13)
  fit2 <- el_sd(women$height, mean = 35, sd = 13, weights = women$weight)
  expect_type(confint(fit, parm = "1"), "double")
  expect_type(confint(fit, parm = c("e", "f")), "double")
  expect_type(confint(fit, parm = c("1", "f")), "double")
  names(fit@coefficients) <- "it"
  expect_type(confint(fit, parm = 1), "double")
  expect_type(confint(fit2, parm = "1"), "double")
})

test_that("`ELMT` objects.", {
  fit <- el_glm(event ~ mag + dist + accel,
                family = poisson("log"),
                data = attenu
  )
  lhs <- list(
    matrix(c(1, 21, 1, 0, 2, 1, 5, 9), nrow = 2),
    matrix(c(0, 1, 0, -1), nrow = 1)
  )
  rhs <- list(c(2, 0), 0)
  out <- elmt(fit, rhs = rhs, lhs = lhs)
  expect_error(confint(out))
  lhs <- list(
    matrix(c(1, 21, 0, 0), nrow = 1),
    matrix(c(0, 1, 0, -1), nrow = 1)
  )
  rhs <- c(2, 0)
  out <- elmt(fit, rhs = rhs, lhs = lhs)
  out@data <- NULL
  expect_error(confint(out))
  out <- elmt(fit, rhs = rhs, lhs = lhs)
  expect_error(confint(out, control = list(maxit = 67)))
  out2 <- confint(out)
  expect_type(out2, "double")
  expect_true(getEstimates(out)[1L] < out2[1L, 2L])
  expect_true(getEstimates(out)[1L] > out2[1L, 1L])
  expect_true(getEstimates(out)[2L] < out2[2L, 2L])
  expect_true(getEstimates(out)[2L] > out2[2L, 1L])
})
