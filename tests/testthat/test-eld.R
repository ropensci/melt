test_that("invalid 'object", {
  skip_on_os("windows", arch = "i386")
  n <- 100
  x <- rnorm(n)
  y <- 1 + x + rnorm(n)
  df <- data.frame(y, x)
  optcfg <- melt_control(maxit = 200L, tol = 1e-08, th = 1e+10)
  fit <- el_lm(y ~ x, df, control = optcfg)
  eld <- eld(fit)
  pdf(NULL)
  plot(eld)
  expect_length(eld, n)
  lhs <- matrix(c(1, -1), nrow = 1)
  fit2 <- lht(fit, lhs = lhs)
  expect_error(eld(df))
  fit$data.matrix <- NULL
  expect_error(eld(fit))
  expect_error(eld(fit2))
})
