# test_that("check convergence", {
#   grid <- seq(-1, 1, length.out = 1000)
#   optcfg <- list(maxit = 20L, abstol = 1e-08, threshold = 1e+10)
#   conv <- function(par) {
#     el_mean(par, x, control = optcfg)$optim$convergence
#   }
#   x <- c(-1.5, 1.5, rnorm(10))
#   expect_identical(all(sapply(grid, conv)), TRUE)
# })
