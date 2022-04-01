test_that("maximum EL estimate", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- (max(x) + min(x)) / 2
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  loglik <- suppressWarnings(logLik(fit, REML = T))
  attributes(loglik) <- NULL
  expect_equal(loglik, -n * log(n))
})

test_that("constrained EL estimate", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.1 * x - x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  lhs <- matrix(c(0, 1, 0), nrow = 1)
  fit2 <- lht(fit, lhs = lhs, control = optcfg)
  expect_output(print(fit2))
  expect_lt(logLik(fit2), -n * log(n))
})
