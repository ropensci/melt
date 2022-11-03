test_that("Invalid `x`.", {
  expect_error(el_sd(1, mean = 0, sd = 1))
  expect_error(el_sd(c(1, Inf), mean = 0, sd = 1))
})

test_that("Invalid `mean`.", {
  x <- sleep$extra
  expect_error(el_sd(x, mean = NA, sd = 1))
  expect_error(el_sd(x, mean = NULL, sd = 1))
  expect_error(el_sd(x, mean = Inf, sd = 1))
  expect_error(el_sd(x, mean = NaN, sd = 1))
  expect_error(el_sd(x, mean = "error", sd = 1))
  expect_error(el_sd(x, mean = c(0, 0), sd = 1))
})

test_that("Invalid `sd`.", {
  x <- sleep$extra
  expect_error(el_sd(x, mean = 0, sd = NA))
  expect_error(el_sd(x, mean = 0, sd = NULL))
  expect_error(el_sd(x, mean = 0, sd = Inf))
  expect_error(el_sd(x, mean = 0, sd = NaN))
  expect_error(el_sd(x, mean = 0, sd = "error"))
  expect_error(el_sd(x, mean = 0, sd = c(0, 0)))
  expect_error(el_sd(x, mean = 0, sd = 0))
  expect_error(el_sd(x, mean = 0, sd = -1))
})

test_that("Probabilities add up to 1.", {
  x <- women$height
  w <- women$weight
  fit <- el_sd(x, mean = 65, sd = 4)
  expect_output(print(fit))
  expect_equal(sum(exp(logProb(fit))), 1, tolerance = 1e-07)
  fit2 <- el_sd(x, mean = 65, sd = 4, weights = w)
  expect_equal(sum(exp(logProb(fit2))), 1, tolerance = 1e-07)
})

test_that("Identical weights means no weights.", {
  x <- women$height
  names(x) <- rownames(women)
  fit <- suppressMessages(el_sd(x,
    mean = 65, sd = 4,
    control = el_control(verbose = TRUE)
  ))
  fit2 <- el_sd(x, mean = 65, sd = 4, weights = rep(0.5, length(x)))
  fit2@weights <- numeric()
  expect_equal(getOptim(fit), getOptim(fit2))
})

test_that("Conversion between `logl` and `loglr`.", {
  x <- women$height
  n <- length(x)
  fit <- el_sd(x, mean = 65, sd = 4)
  expect_equal(logL(fit) + n * log(n), logLR(fit))
  fit2 <- el_sd(x, mean = 65, sd = 4, weights = women$weight)
  w <- weights(fit2)
  expect_equal(logL(fit2) + sum(w * (log(n) - log(w))), logLR(fit2))
})
