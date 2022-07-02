test_that("parallel computation yields the same results (binomial - logit)", {
  skip_on_cran()
  set.seed(54324)
  n <- 100
  p <- 2
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  l <- -1 + x %*% as.vector(b)
  mu <- 1 / (1 + exp(-l))
  y <- vapply(mu, FUN = function(x) rbinom(1, 1, x), FUN.VALUE = integer(1))
  df <- data.frame(cbind(y, x))
  fit <- el_glm(y ~ .,
    family = binomial("logit"), data = df,
    control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = binomial("logit"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
})

test_that("parallel computation yields the same results (binomial - probit)", {
  skip_on_cran()
  set.seed(224543)
  n <- 100
  p <- 2
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  l <- 1.5 + x %*% as.vector(b)
  mu <- pnorm(l)
  y <- vapply(mu, FUN = function(x) rbinom(1, 1, x), FUN.VALUE = integer(1))
  df <- data.frame(cbind(y, x))
  fit <- el_glm(y ~ .,
                family = binomial("probit"), data = df,
                control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = binomial("probit"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
})

test_that("parallel computation yields the same results (binomial - log)", {
  skip_on_cran()
  set.seed(54324)
  n <- 100
  p <- 2
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  l <- -3 + x %*% as.vector(b)
  mu <- exp(l)
  y <- vapply(mu, FUN = function(x) rbinom(1, 1, x), FUN.VALUE = integer(1))
  df <- data.frame(cbind(y, x))
  fit <- el_glm(y ~ .,
    family = binomial("log"), data = df,
    control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = binomial("log"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
})
