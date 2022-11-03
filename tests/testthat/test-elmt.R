test_that("Invalid `rhs`.", {
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  expect_error(elmt(fit))
  expect_error(elmt(fit, rhs = list(c(30, -1, -1))))
  expect_error(elmt(fit, rhs = matrix(c(30, -1, -1), ncol = 3)))
  expect_error(elmt(fit, rhs = Inf))
  expect_error(elmt(fit, rhs = list(c(1, 1, 1), c(1, 1, 1, 1))))
})

test_that("Invalid `lhs`.", {
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  expect_error(elmt(fit, lhs = matrix(c(1, 1, 1), nrow = 1)))
  expect_error(elmt(fit, lhs = list(c(1))))
  expect_error(elmt(fit, lhs = matrix(c(1, 1, 1, 1, 1, NA), nrow = 2)))
  expect_error(elmt(fit, lhs = list(c(1), c(1))))
  expect_error(elmt(fit, lhs = matrix(c(1, 1, 1, 1), nrow = 2)))
  expect_error(elmt(fit, lhs = list(
    matrix(c(1, 2, 3, 4), nrow = 2),
    matrix(c(1, 2, 3, 4), nrow = 2)
  )))
  expect_error(elmt(fit, lhs = matrix(c(0, 1, 0, 1, 0, 1), nrow = 2)))
  expect_error(elmt(fit, lhs = list(
    matrix(c(1, 1, 1, 1, 1, 1), nrow = 2),
    matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  )))
  expect_error(elmt(fit, lhs = list(
    matrix(c(1, 1, 3, 4, 5, 6, 7, 8, 9), nrow = 3),
    matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3)
  )))
})

test_that("Invalid `alpha`.", {
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  rhs <- list(c(0, 0, 0), c(1, 1, 1))
  expect_error(elmt(fit, rhs = rhs, alpha = 1))
})

test_that("Invalid `control`.", {
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  rhs <- c(0, 0, 0)
  lhs <- rbind(c(1, 33, 0), c(0, 1, -100))
  expect_error(elmt(fit, rhs = rhs, lhs = lhs, control = list(maxit = 10)))
})

test_that("Incompatible `rhs` and `lhs`.", {
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  rhs <- c(0, 0, 0)
  lhs <- rbind(c(1, 33, 0), c(0, 1, -100))
  expect_error(elmt(fit, rhs = rhs, lhs = lhs))
})

test_that("`rhs` == `NULL`.", {
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  lhs <- rbind(c(1, 33, 0), c(0, 1, -100))
  out <- elmt(fit, lhs = lhs)
  expect_output(show(out))
  expect_output(print(out))
})

test_that("Matrix `rhs`.", {
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  rhs <- matrix(c(0, 0), ncol = 1)
  lhs <- rbind(c(1, 33, 0), c(0, 1, -100))
  expect_s4_class(elmt(fit, rhs = rhs, lhs = lhs), "ELMT")
})

test_that("List `lhs`.", {
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  lhs <- list(matrix(c(1, 33, 0), nrow = 1), matrix(c(0, 1, -100), nrow = 1))
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
})

test_that("`el_mean()`.", {
  fit <- el_mean(trees, par = c(13, 76, 30))
  lhs <- rbind(c(1, -1, 0), c(0, 1, -1))
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
})

test_that("`el_glm()` (gaussian - inverse).", {
  fit <- el_glm(mpg ~ disp + hp + wt,
    family = gaussian("inverse"),
    data = mtcars, control = el_control(nthreads = 1)
  )
  wfit <- el_glm(mpg ~ disp + hp + wt,
    family = gaussian("inverse"), data = mtcars, weights = gear,
    control = el_control(nthreads = 1)
  )
  lhs <- list(
    matrix(c(0.001, -1, 0, 0), nrow = 1),
    matrix(c(0, 1, -1, 0), nrow = 1)
  )
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
  expect_s4_class(elmt(wfit, lhs = lhs), "ELMT")
})

test_that("`el_glm()` (binomial - logit).", {
  fit <- el_glm(wool ~ ., family = binomial, data = warpbreaks)
  wfit <- el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = breaks
  )
  lhs <- list(
    matrix(c(1, 100, 0, 0), nrow = 1),
    matrix(c(0, 1, -1, 0), nrow = 1)
  )
  out <- elmt(fit, lhs = lhs)
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_output(print(summary(out)))
  expect_output(show(summary(out)))
  expect_equal(sum(exp(logProb(fit))), 1)
  expect_equal(sum(exp(logProb(wfit))), 1)
  lhs2 <- list(
    matrix(c(1, 100, 0, 0, 0, 2, 1, 1), nrow = 2),
    matrix(c(0, 1, -1, 0), nrow = 1)
  )
  out <- elmt(fit, lhs = lhs2)
  expect_output(print(summary(out)))
  expect_output(show(summary(out)))
})

test_that("`el_glm()` (binomial - probit).", {
  fit <- el_glm(wool ~ ., family = binomial("probit"), data = warpbreaks)
  wfit <- el_glm(wool ~ .,
    family = binomial("probit"), data = warpbreaks,
    weights = breaks
  )
  lhs <- list(
    matrix(c(1, 100, 0, 0), nrow = 1),
    matrix(c(0, 1, -1, 0), nrow = 1)
  )
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_equal(sum(exp(fit@logp)), 1)
  expect_equal(sum(exp(wfit@logp)), 1)
})

test_that("`el_glm()` (poisson - log).", {
  fit <- el_glm(event ~ mag + dist + accel,
    family = poisson("log"),
    data = attenu
  )
  wfit <- el_glm(event ~ mag + dist + accel,
    family = poisson("log"),
    data = attenu, weights = mag
  )
  lhs <- list(
    matrix(c(1, 21, 0, 0), nrow = 1),
    matrix(c(0, 1, 0, -1), nrow = 1)
  )
  rhs <- c(2, 0)
  expect_s4_class(elmt(fit, rhs = rhs, lhs = lhs), "ELMT")
  expect_s4_class(elmt(wfit, rhs = rhs, lhs = lhs), "ELMT")
})
