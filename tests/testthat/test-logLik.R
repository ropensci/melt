test_that("maximum EL estimate", {
  skip_on_os("windows", arch = "i386")
  n <- 10
  x <- rnorm(n)
  par <- runif(1, min(x), max(x))
  optcfg <- melt_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  loglik <- suppressWarnings(logLik(fit, REML = T))
  attributes(loglik) <- NULL
  expect_equal(loglik, -n * log(n))
})
