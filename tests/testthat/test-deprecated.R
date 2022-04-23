test_that("el_test", {
  skip_on_os("windows", arch = "i386")
  data("clothianidin")
  expect_warning(out <- el_test(clo ~ trt | blk, clothianidin,
    lhs = matrix(c(
      1, -1, 0, 0,
      0, 1, -1, 0,
      0, 0, 1, -1
    ),
    byrow = TRUE,
    nrow = 3
    )
  ))
  expect_output(print(out))
})
