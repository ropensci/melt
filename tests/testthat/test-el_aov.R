test_that("probabilities add up to 1", {
  skip_on_os("windows", arch = "i386")
  data("clothianidin")
  expect_warning(out <- el_aov(clo ~ trt, clothianidin))
  expect_output(print(out))
})
