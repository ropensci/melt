test_that("Invalid `formula`.", {
  expect_error(el_glm(cbind(group, ID) ~ extra,
    family = binomial,
    data = sleep
  ))
})

test_that("Invalid `family`.", {
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
  expect_error(el_glm(carb ~ .,
    family = quasipoisson("identity"),
    data = mtcars
  ))
})

test_that("Invalid `data`.", {
  x <- warpbreaks$breaks
  x2 <- x
  y <- warpbreaks$wool
  df <- data.frame(y, x, x2)
  expect_error(el_glm(y ~ ., family = binomial, data = df))
})

test_that("Invalid `weights`.", {
  expect_error(el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = rep("error", 54)
  ))
  expect_error(el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = rep(-1, 54)
  ))
})

test_that("Invalid `control`.", {
  expect_error(el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    control = list(maxit = 2L)
  ))
})

test_that("Empty model.", {
  fit <- el_glm(wool ~ 0, family = binomial, data = warpbreaks)
  expect_output(print(summary(fit)))
  fit2 <- el_glm(gear ~ 0, family = quasipoisson("log"), data = mtcars)
  expect_output(print(summary(fit2)))
})

test_that("Probabilities add up to 1.", {
  breaks <- warpbreaks$breaks
  wool <- warpbreaks$wool
  tension <- warpbreaks$tension
  fit <- el_glm(wool ~ breaks + tension, family = binomial)
  wfit <- el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = warpbreaks$breaks
  )
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_equal(sum(exp(fit@logp)), 1)
  expect_equal(sum(exp(wfit@logp)), 1)
})

test_that("conversion between `loglik` and `loglr`.", {
  fit <- el_glm(wool ~ ., family = binomial, data = warpbreaks)
  wfit <- el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = warpbreaks$breaks
  )
  n <- nrow(warpbreaks)
  w <- weights(wfit)
  expect_equal(logL(fit) + n * log(n), logLR(fit))
  expect_equal(logL(wfit) + sum(w * (log(n) - log(w))), logLR(wfit))
})

test_that("No intercept.", {
  fit <- el_glm(wool ~ -1 + ., family = binomial, data = warpbreaks)
  expect_s4_class(fit, "GLM")
})

test_that("`dim` attribute.", {
  wool <- warpbreaks$wool
  dim(wool) <- 54
  breaks <- warpbreaks$breaks
  fit <- el_glm(wool ~ breaks, family = binomial)
  expect_s4_class(fit, "GLM")
})

test_that("`verbose` == TRUE in `el_control()`.", {
  expect_message(el_glm(wool ~ breaks,
    family = binomial, data = warpbreaks, control = el_control(verbose = TRUE)
  ))
})

test_that("`family` == `quasipoisson`.", {
  fit <- el_glm(carb ~ .,
    family = quasipoisson("log"), data = mtcars,
    control = el_control(tol = 1e-04, th = 1000)
  )
  expect_s4_class(fit, "QGLM")
})

test_that("`print()` method.", {
  fit <- el_glm(wool ~ -1 + ., family = binomial, data = warpbreaks)
  out <- summary(fit)
  expect_output(print(out))
  expect_output(show(out))
  out@aliased <- c(TRUE, TRUE, TRUE, TRUE)
  expect_output(print(out))
  df2 <- warpbreaks
  df2[1, 1] <- NA
  fit2 <- el_glm(wool ~ -1 + ., family = binomial, data = df2)
  out2 <- summary(fit2)
  expect_output(print(out2))
  fit3 <- el_glm(gear ~ mpg + cyl + disp,
    family = quasipoisson("log"), data = mtcars
  )
  out3 <- summary(fit3)
  expect_output(print(out3))
  expect_output(show(out3))
  out3@aliased <- c(TRUE, TRUE, TRUE, TRUE)
  expect_output(print(out3))
  df3 <- mtcars
  df3[1, 1] <- NA
  fit4 <- el_glm(gear ~ mpg + cyl + disp,
    family = quasipoisson("log"), data = df3
  )
  out4 <- summary(fit4)
  expect_output(print(out4))
})
