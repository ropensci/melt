test_that("convergence check", {
  skip_on_os("windows", arch = "i386")
  x <- c(-1.5, 1.5, rnorm(10))
  grid <- seq(-1, 1, length.out = 1000)
  conv <- function(par) {
    el_mean(par, x, control = list(maxit = 20L, abstol = 1e-08,
                                    threshold = 1e+10))$optim$convergence
  }
  expect_true(all(sapply(grid, conv)))
})
