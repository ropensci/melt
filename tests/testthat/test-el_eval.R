test_that("check convergence", {
  par <- 1
  x <- c(rnorm(10), 2)
  g <- x - 1
  optcfg <- list(maxit = 10L, abstol = 1e-05, threshold = 1e+01)

  a1 <- el_eval(g, control = optcfg)$optim
  a2 <- el_mean(par, x, control = optcfg)$optim
  a2[["type"]] <- NULL
  expect_identical(a1$convergence, TRUE)
  expect_identical(a1, a2)
})
