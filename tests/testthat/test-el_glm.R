test_that("invalid 'formula'", {
  expect_error(el_glm(cbind(group, ID) ~ extra,
    family = binomial,
    data = sleep
  ))
})

test_that("invalid 'family'", {
  expect_error(el_glm(ID ~ extra, family = binomial("cauchit"), data = sleep))
  expect_error(el_glm(ID ~ extra, family = binomial("cloglog"), data = sleep))
  expect_error(el_glm(ID ~ extra, family = "error", data = sleep))
  expect_error(capture.output(el_glm(y ~ x, family = function() {}, df)))
  expect_error(el_glm(Wind ~ Temp,
    family = gaussian(link = make.link("sqrt")),
    data = airquality
  ))
  expect_error(el_glm(Temp ~ Wind,
    family = poisson(link = make.link("inverse")),
    data = airquality
  ))
  expect_error(el_glm(Wind ~ Temp,
    family = inverse.gaussian("identity"),
    data = airquality
  ))
})

test_that("invalid 'data'", {
  x <- warpbreaks$breaks
  x2 <- x
  y <- warpbreaks$wool
  df <- data.frame(y, x, x2)
  expect_error(el_glm(y ~ ., family = binomial, data = df))
})

test_that("invalid 'weights'", {
  expect_error(el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = rep("error", 54)
  ))
  expect_error(el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = rep(-1, 54)
  ))
})

test_that("invalid 'control'", {
  expect_error(el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    control = list(maxit = 2L)
  ))
})

test_that("empty model", {
  optcfg <- el_control(tol = 1e-08, th = 1e+10)
  fit <- el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    control = optcfg
  )
  expect_output(print(summary(fit)))
})

test_that("probabilities add up to 1", {
  fit <- el_glm(wool ~ ., family = binomial, data = warpbreaks)
  wfit <- el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = warpbreaks$breaks
  )
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_equal(sum(exp(fit@logp)), 1)
  expect_equal(sum(exp(wfit@logp)), 1)
})

test_that("conversion between loglik and loglr", {
  fit <- el_glm(wool ~ ., family = binomial, data = warpbreaks)
  wfit <- el_glm(wool ~ .,
                 family = binomial, data = warpbreaks,
                 weights = warpbreaks$breaks
  )
  n <- nrow(warpbreaks)
  w <- weights(wfit)
  expect_equal(fit@logl + n * log(n), fit@loglr)
  expect_equal(wfit@logl + sum(w * (log(n) - log(w))), wfit@loglr)
})








# test_that("probabilities add up to 1", {
#   skip_on_os("windows", arch = "i386")
#   n <- 50
#   x <- rnorm(n)
#   x2 <- rnorm(n)
#   l <- -2 + 0.2 * x + 1 * x2
#   mu <- 1 / (1 + exp(-l))
#   y <- rbinom(n, 1, mu)
#   df <- data.frame(y, x, x2)
#   optcfg <- el_control(tol = 1e-08, th = 1e+10)
#   fit <- el_glm(y ~ x + x2, family = binomial, df, control = optcfg)
#   expect_output(print(fit))
#   expect_output(print(summary(fit)))
#   expect_equal(sum(exp(fit@logp)), 1)
# })
#
# test_that("probabilities add up to 1 (weighted)", {
#   skip_on_os("windows", arch = "i386")
#   n <- 50
#   x <- rnorm(n)
#   x2 <- rnorm(n)
#   l <- -2 + 0.2 * x + 1 * x2
#   mu <- 1 / (1 + exp(-l))
#   y <- rbinom(n, 1, mu)
#   df <- data.frame(y, x, x2)
#   w <- 1 + runif(n, min = -0.5, max = 0.5)
#   optcfg <- el_control(tol = 1e-08, th = 1e+10)
#   fit <- el_glm(y ~ x + x2, family = binomial, df, control = optcfg)
#   expect_output(print(fit))
#   expect_output(print(summary(fit)))
#   expect_equal(sum(exp(fit@logp)), 1)
# })


# test_that("loglik to loglr", {
#   skip_on_os("windows", arch = "i386")
#   n <- 50
#   x <- rnorm(n)
#   x2 <- rnorm(n)
#   l <- -1 + 0.9 * x + 0.3 * x2
#   mu <- 1 / (1 + exp(-l))
#   y <- rbinom(n, 1, mu)
#   df <- data.frame(y, x, x2)
#   optcfg <- el_control(tol = 1e-08, th = 1e+10)
#   fit <- el_glm(y ~ x + x2, family = binomial, df, control = optcfg)
#   expect_equal(fit@logl + n * log(n), fit@loglr)
# })
#
# test_that("loglik to loglr (weighted)", {
#   skip_on_os("windows", arch = "i386")
#   n <- 50
#   x <- rnorm(n)
#   x2 <- rnorm(n)
#   l <- -1 + 0.9 * x + 0.3 * x2
#   mu <- 1 / (1 + exp(-l))
#   y <- rbinom(n, 1, mu)
#   df <- data.frame(y, x, x2)
#   w <- 1 + runif(n, min = -0.5, max = 0.5)
#   optcfg <- el_control(tol = 1e-08, th = 1e+10)
#   fit <- suppressWarnings(el_glm(y ~ x + x2,
#     family = binomial, df,
#     weights = w, control = optcfg
#   ))
#   w <- fit@weights
#   expect_equal(fit@logl + sum(w * (log(n) - log(w))), fit@loglr)
# })
