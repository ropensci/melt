test_that("probabilities add up to 1", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  expect_equal(sum(exp(fit$log.prob)), 1)
})

# test_that("probabilities add up to 1 (weighted)", {
#   skip_on_os("windows", arch = "i386")
#   n <- 20
#   x <- rnorm(n)
#   x2 <- rnorm(n)
#   y <- 1 + 0.1 * x - 0.1 * x2 + rnorm(n)
#   df <- data.frame(y, x, x2)
#   w <- 1 + runif(20, min = -0.01, max = 0.01)
#   optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
#   fit <- el_lm(y ~ x + x2, df, weights =  w, control = optcfg)
#   expect_equal(sum(exp(fit$log.prob)), 1)
# })
