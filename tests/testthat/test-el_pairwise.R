test_that("output", {
  skip_on_cran()
  data("clothianidin")
  out1 <- el_pairwise(clo ~ trt | blk, clothianidin, B = 500, progress = TRUE)
  out2 <- el_pairwise(clo ~ trt | blk, clothianidin,
    control = "Naked",
    method = "NB", progress = TRUE, B = 2000, nthread = 2
  )
  expect_output(print(out1))
  expect_output(print(out2))
})
