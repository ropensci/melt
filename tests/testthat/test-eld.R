test_that("Invalid `object`.", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  fit@data <- matrix(NA_real_, nrow = 0L, ncol = 0L)
  expect_error(eld(fit))
})

test_that("Invalid `control`.", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(eld(fit, control = list(maxit = 20L)))
})

test_that("`weights`", {
  fit <- el_lm(extra ~ ID, data = sleep, weights = as.numeric(group))
  expect_s4_class(eld(fit), "ELD")
})

test_that("Mean parameter.", {
  fit <- el_mean(women$height, par = 67)
  fit2 <- el_mean(women$height, par = 67, weights = women$weight)
  expect_s4_class(eld(fit), "ELD")
  expect_s4_class(eld(fit2), "ELD")
})

test_that("`SD` class.", {
  x <- women$height
  w <- women$weight
  fit <- el_sd(x, mean = 65, sd = 4)
  wfit <- el_sd(x, mean = 65, sd = 4, weights = w)
  expect_s4_class(eld(fit), "ELD")
  expect_s4_class(eld(wfit), "ELD")
})

test_that("`GLM` class.", {
  fit <- el_glm(wool ~ -1 + breaks, family = binomial, data = warpbreaks)
  wfit <- el_glm(wool ~ .,
    family = binomial, data = warpbreaks,
    weights = breaks
  )
  expect_s4_class(eld(fit), "ELD")
  expect_s4_class(eld(wfit), "ELD")
})

test_that("`plot()` method.", {
  fit <- el_lm(mpg ~ disp + hp + wt, data = mtcars)
  eld <- eld(fit)
  pdf(NULL)
  expect_invisible(plot(fit))
  expect_invisible(plot(eld))
})

test_that("Missing `object`.", {
  expect_null(eld())
})
