test_that("deprecated", {
  n <- 10
  x <- rnorm(n)
  par <- runif(1, min(x), max(x))
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_mean(par, x, control = optcfg)
  expect_warning(lht(fit, lhs = 1, control = optcfg))
})

test_that("deprecated", {
  skip_on_cran()
  data("clothianidin")
  expect_warning(
    out1 <- el_pairwise(clo ~ trt | blk, clothianidin, B = 500, progress = TRUE)
  )
  expect_warning(
    out2 <- el_pairwise(clo ~ trt | blk, clothianidin,
      control = "Naked",
      method = "NB", progress = TRUE, B = 2000, nthread = 2
    )
  )
  expect_output(print(out1))
  expect_output(print(out2))
})
