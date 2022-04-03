test_that("invalid 'object", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- (max(x) + min(x)) / 2
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, w, control = optcfg)
  fit$data.matrix <- NULL
  expect_error(lht(fit, lhs = 1))
  class(fit) <- NULL
  expect_error(lht(fit, lhs = 1))
})

test_that("invalid 'lhs' and 'rhs'", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x + x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  lhs <- matrix(c(0, 1, -1), nrow = 1)
  expect_error(lht(fit, c(NA, 0), lhs))
  expect_error(lht(fit, c(1, 0), lhs))

  w <- 1 + runif(n, min = -0.5, max = 0.5)
  fit2 <- el_lm(y ~ x + x2, df, weights =  w, control = optcfg)
  out <- lht(fit2, lhs = lhs, control = optcfg)
  expect_lt(out$statistic, 20)

  lhs3 <- matrix(c(1, 1, 0, 0, 0, 0), nrow = 2)
  expect_error(lht(fit, lhs = lhs3))
  expect_error(lht(fit2, lhs = lhs3))
})
