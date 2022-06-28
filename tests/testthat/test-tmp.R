test_that("same results with parallel computing", {
  fit <- el_lm(mpg ~ ., data = mtcars, control = el_control(nthreads = 1))
  fit2 <- el_lm(mpg ~ ., data = mtcars)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
  wfit <- el_lm(mpg ~ .,
                data = mtcars, weights = mtcars$wt,
                control = el_control(nthreads = 1)
  )
  wfit2 <- el_lm(mpg ~ ., data = mtcars, weights = mtcars$wt)
  expect_equal(wfit@optim, wfit2@optim)
  expect_equal(wfit@parTests, wfit2@parTests)
})

test_that("same results with parallel computing (binomial)", {
  skip_on_os("windows", arch = "i386")
  # skip_on_ci()
  skip_on_cran()
  n <- 500
  p <- 15
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  l <- 1 + x %*% as.vector(b)
  mu <- 1 / (1 + exp(-l))
  y <- rbinom(n, 1, mu)
  w <- 1 + runif(n, min = -0.5, max = 0.5)
  df <- data.frame(cbind(y, x))
  lfit <- el_glm(y ~ .,
                 family = binomial(link = "logit"), df,
                 control = el_control(tol = 1e-08, th = 1e+10, nthreads = 1)
  )
  lfit2 <- el_glm(y ~ .,
                  family = binomial(link = "logit"), df,
                  control = el_control(tol = 1e-08, th = 1e+10)
  )
  expect_equal(lfit@optim, lfit2@optim)
  expect_equal(lfit@parTests, lfit2@parTests)
  wlfit <- suppressWarnings(
    el_glm(y ~ .,
           family = binomial(link = "logit"), df, weights = w,
           control = el_control(tol = 1e-08, th = 1e+10, nthreads = 1)
    )
  )
  wlfit2 <- suppressWarnings(
    el_glm(y ~ .,
           family = binomial(link = "logit"), df, weights = w,
           control = el_control(tol = 1e-08, th = 1e+10)
    )
  )
  expect_equal(wlfit@optim, wlfit2@optim)
  expect_equal(wlfit@parTests, wlfit2@parTests)

  pfit <- el_glm(y ~ .,
                 family = binomial(link = "probit"), df,
                 control = el_control(tol = 1e-08, th = 1e+10, nthreads = 1)
  )
  pfit2 <- el_glm(y ~ .,
                  family = binomial(link = "probit"), df,
                  control = el_control(tol = 1e-08, th = 1e+10)
  )
  expect_equal(pfit@optim, pfit2@optim)
  expect_equal(pfit@parTests, pfit2@parTests)
  wpfit <- suppressWarnings(
    el_glm(y ~ .,
           family = binomial(link = "probit"), df, weights = w,
           control = el_control(tol = 1e-08, th = 1e+10, nthreads = 1)
    )
  )
  wpfit2 <- suppressWarnings(
    el_glm(y ~ .,
           family = binomial(link = "probit"), df, weights = w,
           control = el_control(tol = 1e-08, th = 1e+10)
    )
  )
  expect_equal(wpfit@optim, wpfit2@optim)
  expect_equal(wpfit@parTests, wpfit2@parTests)
})
