data("mtcars")
fit <- el_lm(mpg ~ cyl + disp, mtcars)

test_that("invalid 'rhs'", {
  # no hypothesis
  expect_error(elmt(fit))
  # single hypothesis
  expect_error(elmt(fit, rhs = list(c(30, -1, -1))))
  # matrix rhs
  expect_error(elmt(fit, rhs = matrix(c(30, -1, -1), ncol = 3)))
  # non-numeric rhs
  expect_error(elmt(fit, rhs = Inf))
  # length
  expect_error(elmt(fit, rhs = list(c(1, 1, 1), c(1, 1, 1, 1))))
})

test_that("invalid 'lhs'", {
  # single hypothesis
  expect_error(elmt(fit, lhs = matrix(c(1, 1, 1), nrow = 1)))
  expect_error(elmt(fit, lhs = list(c(1))))
  # non-numeric matrix
  expect_error(elmt(fit, lhs = matrix(c(1, 1, 1, 1, 1, NA), nrow = 2)))
  expect_error(elmt(fit, lhs = list(c(1), c(1))))
  # length
  expect_error(elmt(fit, lhs = matrix(c(1, 1, 1, 1), nrow = 2)))
  expect_error(elmt(fit, lhs = list(
    matrix(c(1, 2, 3, 4), nrow = 2),
    matrix(c(1, 2, 3, 4), nrow = 2)
  )))
  # zero vector
  expect_error(elmt(fit, lhs = matrix(c(0, 1, 0, 1, 0, 1), nrow = 2)))
  # rank
  expect_error(elmt(fit, lhs = list(
    matrix(c(1, 1, 1, 1, 1, 1), nrow = 2),
    matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  )))
  expect_error(elmt(fit, lhs = list(
    matrix(c(1, 1, 3, 4, 5, 6, 7, 8, 9), nrow = 3),
    matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3)
  )))
})

test_that("invalid 'alpha'", {
  rhs <- list(c(0, 0, 0), c(1, 1, 1))
  expect_error(elmt(fit, rhs = rhs, alpha = 1))
})

test_that("'rhs' == NULL", {
  lhs <- rbind(c(1, 33, 0), c(0, 1, -100))
  expect_s4_class(elmt(fit, lhs = lhs), "MELT")
})
