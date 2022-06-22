test_that("deprecated", {
  n <- 10
  x <- rnorm(n)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_mean(x, par, control = optcfg)
  expect_warning(lht(fit, lhs = 1, control = optcfg))
})
