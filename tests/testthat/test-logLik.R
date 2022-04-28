test_that("logLik at maximum EL estimates", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  loglik <- suppressWarnings(logLik(fit, REML = T))
  expect_equal(loglik@logLik, -n * log(n))
})

test_that("empty model", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  y <- 1 + x + rnorm(n)
  df <- data.frame(x, y)
  fit <- el_lm(y ~ 0, df)
  expect_error(logLik(fit))
})
