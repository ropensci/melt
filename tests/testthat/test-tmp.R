# test_that("same results with parallel computing", {
#   fit <- el_lm(mpg ~ ., data = mtcars, control = el_control(nthreads = 1))
#   fit2 <- el_lm(mpg ~ ., data = mtcars)
#   expect_equal(fit@optim, fit2@optim)
#   expect_equal(fit@parTests, fit2@parTests)
#   wfit <- el_lm(mpg ~ .,
#                 data = mtcars, weights = mtcars$wt,
#                 control = el_control(nthreads = 1)
#   )
#   wfit2 <- el_lm(mpg ~ ., data = mtcars, weights = mtcars$wt)
#   expect_equal(wfit@optim, wfit2@optim)
#   expect_equal(wfit@parTests, wfit2@parTests)
# })
