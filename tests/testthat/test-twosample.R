test_that("two sample test", {
  set.seed(1)
  x <- rnorm(100)
  y <- rnorm(100)
  expect_equal(test2sample2_cpp(x,y)$convergence, 1)
})
