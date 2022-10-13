test_that("Invalid `formula`.", {
  expect_error(el_lm(cbind(mtcars$mpg, mtcars$cyl) ~ 1, data = mtcars))
})

test_that("Invalid `data`.", {
  x <- mtcars$disp
  x2 <- mtcars$cyl
  y <- mtcars$mpg
  df <- data.frame(y, x, x2)
  error <- function(x) {}
  df <- list(y, x, x2, error)
  df2 <- list(y, x, x2, matrix(rnorm(100)))
  df3 <- list(y, x, x2, list())
  expect_error(el_lm(y ~ x + x2, data = df))
  expect_error(el_lm(y ~ x + x2, data = df2))
  expect_error(el_lm(y ~ x + x2, data = df3))
})

test_that("Invalid `offset`.", {
  expect_error(el_lm(height ~ weight, data = women, offset = as.matrix(women)))
  fit <- el_lm(height ~ weight, data = women, offset = weight)
  expect_s4_class(fit, "LM")
})

test_that("Invalid `weights`.", {
  w <- women$weight
  w[1] <- -1
  expect_error(el_lm(height ~ weight, data = women, weights = w))
  w[1] <- "error"
  expect_error(el_lm(height ~ weight, data = women, weights = w))
})

test_that("Invalid `control`.", {
  expect_error(el_lm(mpg ~ ., data = mtcars, control = list()))
})

test_that("Probabilities add up to 1.", {
  x <- women$weight
  y <- 10 + 0.001 * x + rep(c(1, 0.2, 0.5, 2, -1.2), times = 3)
  df <- data.frame(x, y)
  fit <- el_lm(y ~ x, data = df)
  expect_equal(sum(exp(logProb(fit))), 1)
  fit2 <- el_lm(y ~ x, data = df, weights = women$height)
  expect_equal(sum(exp(logProb(fit2))), 1)
})

test_that("Conversion between `logl` and `loglr`.", {
  fit <- el_lm(eruptions ~ waiting, data = faithful)
  n <- nrow(faithful)
  expect_equal(logL(fit) + n * log(n), logLR(fit))
  wfit <- el_lm(eruptions ~ waiting,
    data = faithful,
    weights = faithful$waiting
  )
  w <- weights(wfit)
  expect_equal(logL(wfit) + sum(w * (log(n) - log(w))), logLR(wfit))
})

test_that("Empty model.", {
  fit <- el_lm(eruptions ~ 0, faithful)
  expect_output(print(summary(fit)))
})

test_that("No intercept.", {
  fit <- el_lm(mpg ~ -1 + ., data = mtcars)
  expect_s4_class(fit, "LM")
})

test_that("`print()` method.", {
  fit <- el_lm(mpg ~ disp + hp, data = mtcars)
  expect_output(show(fit))
  expect_output(print(fit))
  out <- summary(fit)
  expect_output(print(out))
  expect_output(show(out))
  out@aliased <- c(TRUE, TRUE, TRUE)
  expect_output(print(out))
  fit@statistic <- numeric()
  expect_output(print(fit))
  df2 <- mtcars
  df2[1, 1] <- NA
  fit2 <- el_lm(mpg ~ disp + hp, data = df2)
  expect_output(print(summary(fit2)))
  fit3 <- el_lm(mpg ~ -1 + disp + hp, data = mtcars)
  expect_output(print(fit3))
  expect_output(print(summary(fit3)))
  wfit <- el_lm(mpg ~ disp + hp, data = mtcars, weights = gear)
  expect_output(print(wfit))
  wout <- summary(wfit)
  expect_output(print(wout))
})

test_that("`verbose` == TRUE in `el_control()`.", {
  expect_message(el_lm(mpg ~ disp + hp,
    data = mtcars, control = el_control(verbose = TRUE)
  ))
})

test_that("`conv()`, `formula()`, and `nobs()` methods.", {
  fit <- el_lm(mpg ~ disp + hp,
    data = mtcars,
    control = el_control(th = 1e+10)
  )
  fit2 <- el_lm(mpg ~ disp + hp,
    data = mtcars,
    control = el_control(step = .Machine$double.eps, th = 1e+10)
  )
  expect_true(conv(fit))
  expect_false(conv(fit2))
  expect_type(formula(fit), "language")
  expect_identical(nobs(fit), nrow(mtcars))
})

test_that("Exact relationship between predictors.", {
  x <- mtcars$disp
  x2 <- x
  y <- mtcars$mpg
  df <- data.frame(y, x, x2)
  w <- mtcars$carb
  expect_error(el_lm(y ~ x + x2, data = df))
  expect_error(el_lm(y ~ x + x2, data = df, weights = w))
})

test_that("All row and column names are preserved.", {
  wfit <- el_lm(mpg ~ -1 + disp + hp, data = mtcars, weights = qsec)
  row_names <- rownames(mtcars)
  column_names <- colnames(mtcars)
  expect_identical(names(wfit@logp), row_names)
  expect_identical(rownames(getData(wfit)), row_names)
  expect_identical(names(weights(wfit)), row_names)
  expect_identical(formula(wfit), formula(wfit@terms))
  expect_identical(nobs(wfit), nrow(mtcars))
  expect_true(all(names(sigTests(wfit)$statistic) %in% column_names))
  expect_true(all(names(sigTests(wfit)$convergence) %in% column_names))
  expect_true(all(names(getOptim(wfit)$par) %in% column_names))
  expect_true(all(colnames(getData(wfit)[, -c(1L, 2L)]) %in% column_names))
  expect_true(all(names(coef(wfit)) %in% column_names))
  expect_true(all(rownames(confint(wfit)) %in% column_names))
})
