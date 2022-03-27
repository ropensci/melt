test_that("convergence check", {
  skip_on_os("windows", arch = "i386")
  x <- c(-1.5, 1.5, rnorm(10))
  grid <- seq(-1, 1, length.out = 1000)
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  conv <- function(par) {
    el_mean(par, x, control = optcfg)$optim$convergence
  }
  expect_true(all(sapply(grid, conv)))
})

test_that("identical weights == no weights", {
  skip_on_os("windows", arch = "i386")
  x <- rnorm(10)
  par <- runif(1, min(x), max(x))
  w <- rep(runif(1), length(x))
  optcfg <- list(maxit = 200L, tol = 1e-08, th = 1e+10)
  a1 <- el_mean(par, x, control = optcfg)$optim
  a2 <- el_mean(par, x, w, control = optcfg)$optim
  expect_equal(a1$lambda, a2$lambda)
  expect_equal(a1$logLR, a2$logLR)
})
