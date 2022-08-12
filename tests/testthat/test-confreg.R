test_that("Invalid `object`.", {
  fit <- el_lm(mpg ~ 0, data = mtcars)
  fit2 <- el_glm(gear ~ 0, family = quasipoisson("log"), data = mtcars)
  fit3 <- el_glm(gear ~ 1, family = quasipoisson("log"), data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2)))
  expect_error(confreg(fit2, parm = c(1, 2)))
  expect_error(confreg(fit3, parm = c(1, 2)))
  fit@data <- matrix()
  expect_error(confreg(fit, parm = c(1, 2)))
  fit2@data <- matrix()
  expect_error(confreg(fit2, parm = c(1, 2)))
})

test_that("Invalid `parm`.", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2, 3)))
  expect_error(confreg(fit, parm = c(1, 1)))
  expect_error(confreg(fit, parm = c("error", "error2")))
  expect_error(confreg(fit, parm = c(NaN, NA)))
  names(fit@coefficients) <- NULL
  expect_error(confreg(fit, parm = c("error", "error2")))
  fit2 <- el_mean(women$height, 67)
  expect_error(confreg(fit2, parm = c(1, 2)))
  fit2 <- el_glm(gear ~ mpg, family = quasipoisson("log"), data = mtcars)
  expect_error(confreg(fit2, parm = c(1, 2, 3)))
  expect_error(confreg(fit2, parm = c(1, 1)))
  expect_error(confreg(fit2, parm = c("error", "error2")))
  expect_error(confreg(fit2, parm = c(NaN, NA)))
  names(fit2@coefficients) <- NULL
  expect_error(confreg(fit2, parm = c("error", "error2")))
})

test_that("Invalid `level`.", {
  fit <- el_lm(mpg ~ hp + disp, data = mtcars)
  expect_s4_class(confreg(fit, parm = c(1, 2), level = 0), "ConfregEL")
  expect_error(confreg(fit, parm = c(1, 2), level = 1))
  fit2 <- el_glm(gear ~ mpg + cyl, family = quasipoisson("log"), data = mtcars)
  expect_s4_class(confreg(fit2, parm = c(1, 2), level = 0), "ConfregEL")
  expect_error(confreg(fit2, parm = c(1, 2), level = 1))
})

test_that("Invalid `cv`.", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2), cv = 1e+10))
  expect_error(confreg(fit,
    parm = c(1, 2), cv = 1e+10, control = el_control(th = 100)
  ))
  fit2 <- el_glm(gear ~ mpg, family = quasipoisson("log"), data = mtcars)
  expect_error(confreg(fit2, parm = c(1, 2), cv = 1e+10))
  expect_error(confreg(fit2,
    parm = c(1, 2), cv = 1e+10, control = el_control(th = 100)
  ))
})

test_that("Invalid `npoints`.", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2), npoints = -10))
  fit2 <- el_glm(gear ~ 0, family = quasipoisson("log"), data = mtcars)
  expect_error(confreg(fit2, parm = c(1, 2), npoints = -10))
})

test_that("Invalid `control`.", {
  fit <- el_lm(mpg ~ hp, data = mtcars)
  expect_error(confreg(fit, parm = c(1, 2), control = list(maxit = 10)))
  fit2 <- el_glm(gear ~ 0, family = quasipoisson("log"), data = mtcars)
  expect_error(confreg(fit2, parm = c(1, 2), control = list()))
})

test_that("`plot()` method.", {
  fit <- el_lm(mpg ~ disp + hp, data = mtcars)
  cr <- confreg(fit)
  cr2 <- confreg(fit, parm = c(1, 2), cv = 6)
  pdf(NULL)
  expect_invisible(plot(cr2))
  fit2 <- el_glm(gear ~ mpg + cyl,
    family = quasipoisson("log"), data = mtcars
  )
  cr3 <- confreg(fit2)
  cr4 <- confreg(fit2, parm = c(1, 2), cv = 6)
  pdf(NULL)
  expect_invisible(plot(cr4))
})
