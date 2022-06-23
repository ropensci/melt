test_that("invalid 'alpha'", {
  data("mtcars")
  optcfg <- el_control(maxit_l = 200L, tol_l = 1e-08, th = 1e+10)
  fit <- el_lm(mpg ~ cyl + disp, mtcars, control = optcfg)
  rhs <- list(c(0, 0, 0), c(1, 1, 1))
  expect_error(elmt(fit, rhs = rhs, alpha = 1))
})
