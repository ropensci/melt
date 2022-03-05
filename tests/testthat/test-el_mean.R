test_that("check convergence", {
  par <- 1
  x <- c(rnorm(10), 2)
  optcfg <- list(maxit = 10L, abstol = 1e-05, threshold = 1e+01)
  expect_equal(el_mean(par, x, control = optcfg)$optim$convergence, TRUE)
})
