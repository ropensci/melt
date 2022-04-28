test_that("normal", {
  n <- 100
  x <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.4 * x - 0.2 * x2 + rnorm(n)
  df <- data.frame(y, x, x2)
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x + x2, df, control = optcfg)
  cr <- confreg(fit, parm = c(1, 2))
  pdf(NULL)
  plot(cr)
  expect_true(is(cr, "ConfregEL"))
})
