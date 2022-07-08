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
  expect_equal(sum(exp(fit@logp)), 1)
  fit2 <- el_lm(y ~ x, data = df, weights = women$height)
  expect_equal(sum(exp(fit2@logp)), 1)
})

test_that("Conversion between `loglik` and `loglr`.", {
  fit <- el_lm(eruptions ~ waiting, data = faithful)
  n <- nrow(faithful)
  expect_equal(fit@logl + n * log(n), logLR(fit))
  wfit <- el_lm(eruptions ~ waiting,
    data = faithful,
    weights = faithful$waiting
  )
  w <- weights(wfit)
  expect_equal(wfit@logl + sum(w * (log(n) - log(w))), logLR(wfit))
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
})

test_that("`verbose` == TRUE in `el_control()`.", {
  expect_message(el_lm(mpg ~ disp + hp,
    data = mtcars, control = el_control(verbose = TRUE)
  ))
})

test_that("`conv()`, `formula()`, and `nobs()` methods.", {
  fit <- el_lm(mpg ~ disp + hp, data = mtcars)
  expect_false(conv(fit))
  expect_type(formula(fit), "language")
  expect_identical(nobs(fit), nrow(mtcars))
})

#' @srrstats {G5.6, G5.6a, G5.6b} `el_lm()` returns the expected coefficients
#'   (`rep(1, p)`) for a simulated data set generated from a linear model. The
#'   parameters used are `n = 1e+05`, `p = 3`, and `tolerance = 1e-02`, with
#'   three different seeds.
test_that("Parameter recovery tests.", {
  set.seed(5524325)
  n <- 1e+05
  p <- 3
  e <- rnorm(n)
  b <- rep(1, p)
  x <- matrix(rnorm(n * p), ncol = p)
  y <- x %*% b + e
  df <- data.frame(y, x)
  fit <- el_lm(y ~ -1 + ., df)
  expect_equal(max(abs(coef(fit))), 1, tolerance = 1e-02)
  set.seed(55267654)
  x2 <- matrix(rnorm(n * p), ncol = p)
  y2 <- x2 %*% b + e
  df2 <- data.frame(y2, x2)
  fit2 <- el_lm(y2 ~ -1 + ., df2)
  expect_equal(max(abs(coef(fit2))), 1, tolerance = 1e-02)
  set.seed(0841)
  x3 <- matrix(rnorm(n * p), ncol = p)
  y3 <- x3 %*% b + e
  df3 <- data.frame(y3, x3)
  fit3 <- el_lm(y3 ~ -1 + ., df3)
  expect_equal(max(abs(coef(fit3))), 1, tolerance = 1e-02)
})

#' @srrstats {RE7.0, RE7.0a} Exact relationship between predictors causes erros.
test_that("Exact relationship between predictors.", {
  x <- mtcars$disp
  x2 <- x
  y <- mtcars$mpg
  df <- data.frame(y, x, x2)
  w <- mtcars$carb
  expect_error(el_lm(y ~ x + x2, data = df))
  expect_error(el_lm(y ~ x + x2, data = df, weights = w))
})

#' @srrstats {RE7.2, RE7.3} All row or column names (if any) are preserved in
#'   the fitted objects and can be retrieved by accessor methods.
test_that("All row and column names are preserved.", {
  wfit <- el_lm(mpg ~ -1 + disp + hp, data = mtcars, weights = qsec)
  row_names <- rownames(mtcars)
  column_names <- colnames(mtcars)
  expect_identical(names(wfit@logp), row_names)
  expect_identical(rownames(wfit@data), row_names)
  expect_identical(names(weights(wfit)), row_names)
  expect_identical(formula(wfit), formula(wfit@terms))
  expect_identical(nobs(wfit), nrow(mtcars))
  expect_true(all(names(wfit@parTests$statistic) %in% column_names))
  expect_true(all(names(wfit@parTests$convergence) %in% column_names))
  expect_true(all(names(wfit@optim$par) %in% column_names))
  expect_true(all(colnames(wfit@data[, -1L]) %in% column_names))
  expect_true(all(names(coef(wfit)) %in% column_names))
  expect_true(all(rownames(confint(wfit)) %in% column_names))
})
