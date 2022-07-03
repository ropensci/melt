fit <- el_lm(mpg ~ cyl + disp, data = mtcars)

test_that("Invalid `rhs`.", {
  expect_error(elmt(fit))
  expect_error(elmt(fit, rhs = list(c(30, -1, -1))))
  expect_error(elmt(fit, rhs = matrix(c(30, -1, -1), ncol = 3)))
  expect_error(elmt(fit, rhs = Inf))
  expect_error(elmt(fit, rhs = list(c(1, 1, 1), c(1, 1, 1, 1))))
})

test_that("Invalid `lhs`.", {
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
  rhs <- list(c(0, 0, 0), c(1, 1, 1))
  expect_error(elmt(fit, rhs = rhs, alpha = 1))
})

test_that("Incompatible `rhs` and `lhs`.", {
  rhs <- c(0, 0, 0)
  lhs <- rbind(c(1, 33, 0), c(0, 1, -100))
  expect_error(elmt(fit, rhs = rhs, lhs = lhs))
})

test_that("`rhs` == `NULL`.", {
  lhs <- rbind(c(1, 33, 0), c(0, 1, -100))
  out <- elmt(fit, lhs = lhs)
  expect_output(show(out))
  expect_output(print(out))
})

test_that("Matrix `rhs`.", {
  rhs <- matrix(c(0, 0), ncol = 1)
  lhs <- rbind(c(1, 33, 0), c(0, 1, -100))
  expect_s4_class(elmt(fit, rhs = rhs, lhs = lhs), "ELMT")
})

test_that("List `lhs`.", {
  lhs <- list(matrix(c(1, 33, 0), nrow = 1), matrix(c(0, 1, -100), nrow = 1))
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
})
