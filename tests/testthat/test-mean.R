test_that("test mean", {
  set.seed(1)
  x <- matrix(rnorm(100))
  expect_equal(el_mean2(0, x)$convergence, TRUE)
})
