test_that("Parallel computation yields the same results (gaussian - log).", {
  skip_on_cran()
  set.seed(23)
  n <- 100
  p <- 3
  b <- rnorm(p, mean = 0, sd = 0.1)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 2), times = 50)
  l <- 3 + x %*% as.vector(b)
  mu <- exp(l)
  y <- vapply(mu, FUN = function(x) rnorm(1, x), FUN.VALUE = numeric(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ .,
    family = gaussian("log"), data = df,
    control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = gaussian("log"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
  wfit <- el_glm(y ~ .,
    family = gaussian("log"), data = df, weights = w,
    control = el_control(nthreads = 1)
  )
  wfit2 <- el_glm(y ~ ., family = gaussian("log"), data = df, weights = w)
  expect_equal(wfit@optim, wfit2@optim)
  expect_equal(wfit@parTests, wfit2@parTests)
  lhs <- list(
    matrix(c(0, 0, 2, 1.2), nrow = 1),
    matrix(c(0, 1, 0, -.5), nrow = 1)
  )
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
})

test_that("Parallel computation yields the same results (binomial - logit).", {
  skip_on_cran()
  set.seed(54324)
  n <- 100
  p <- 2
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 2), times = 50)
  l <- -1 + x %*% as.vector(b)
  mu <- 1 / (1 + exp(-l))
  y <- vapply(mu, FUN = function(x) rbinom(1, 1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ .,
    family = binomial("logit"), data = df,
    control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = binomial("logit"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
  wfit <- el_glm(y ~ .,
    family = binomial("logit"), data = df, weights = w,
    control = el_control(nthreads = 1)
  )
  wfit2 <- el_glm(y ~ ., family = binomial("logit"), data = df, weights = w)
  expect_equal(wfit@optim, wfit2@optim)
  expect_equal(wfit@parTests, wfit2@parTests)
})

test_that("Parallel computation yields the same results (binomial - probit).", {
  skip_on_cran()
  set.seed(224543)
  n <- 100
  p <- 2
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 2), times = 50)
  l <- 1.5 + x %*% as.vector(b)
  mu <- pnorm(l)
  y <- vapply(mu, FUN = function(x) rbinom(1, 1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ .,
    family = binomial("probit"), data = df,
    control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = binomial("probit"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
  wfit <- el_glm(y ~ .,
    family = binomial("probit"), data = df, weights = w,
    control = el_control(nthreads = 1)
  )
  wfit2 <- el_glm(y ~ ., family = binomial("probit"), data = df, weights = w)
  expect_equal(wfit@optim, wfit2@optim)
  expect_equal(wfit@parTests, wfit2@parTests)
})

test_that("Parallel computation yields the same results (binomial - log).", {
  skip_on_cran()
  set.seed(54324)
  n <- 100
  p <- 2
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 1.0005), times = 50)
  l <- -3 + x %*% as.vector(b)
  mu <- exp(l)
  y <- vapply(mu, FUN = function(x) rbinom(1, 1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ .,
    family = binomial("log"), data = df,
    control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = binomial("log"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
  wfit <- el_glm(y ~ .,
    family = binomial("log"), data = df, weights = w,
    control = el_control(nthreads = 1)
  )
  wfit2 <- el_glm(y ~ ., family = binomial("log"), data = df, weights = w)
  expect_equal(wfit@optim, wfit2@optim)
  expect_equal(wfit@parTests, wfit2@parTests)
})

test_that("Parallel computation yields the same results (poisson - log).", {
  skip_on_cran()
  set.seed(525)
  n <- 100
  p <- 3
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 2), times = 50)
  l <- -2 + x %*% as.vector(b)
  mu <- exp(l)
  y <- vapply(mu, FUN = function(x) rpois(1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ .,
    family = poisson("log"), data = df,
    control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = poisson("log"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
  wfit <- el_glm(y ~ .,
    family = poisson("log"), data = df, weights = w,
    control = el_control(nthreads = 1)
  )
  wfit2 <- el_glm(y ~ ., family = poisson("log"), data = df, weights = w)
  expect_equal(wfit@optim, wfit2@optim)
  expect_equal(wfit@parTests, wfit2@parTests)
})

test_that(
  "Parallel computation yields the same results (poisson - identity).",
  {
    skip_on_cran()
    set.seed(534)
    n <- 100
    p <- 3
    b <- rnorm(p, sd = 0.5)
    x <- matrix(rnorm(n * p), ncol = p)
    w <- rep(c(1, 2), times = 50)
    l <- 3 + x %*% as.vector(b)
    mu <- l
    y <- vapply(mu, FUN = function(x) rpois(1, x), FUN.VALUE = integer(1))
    df <- data.frame(y, x)
    fit <- el_glm(y ~ .,
      family = poisson("identity"), data = df,
      control = el_control(nthreads = 1)
    )
    fit2 <- el_glm(y ~ ., family = poisson("identity"), data = df)
    expect_equal(fit@optim, fit2@optim)
    expect_equal(fit@parTests, fit2@parTests)
    wfit <- el_glm(y ~ .,
      family = poisson("identity"), data = df, weights = w,
      control = el_control(nthreads = 1)
    )
    wfit2 <- el_glm(y ~ ., family = poisson("identity"), data = df, weights = w)
    expect_equal(wfit@optim, wfit2@optim)
    expect_equal(wfit@parTests, wfit2@parTests)
  }
)

test_that("Parallel computation yields the same results (poisson - sqrt).", {
  skip_on_cran()
  set.seed(15234)
  n <- 100
  p <- 3
  b <- rnorm(p, mean = 0, sd = .2)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 2), times = 50)
  l <- 0.5 + x %*% as.vector(b)
  mu <- l^2
  y <- vapply(mu, FUN = function(x) rpois(1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ .,
    family = poisson("sqrt"), data = df,
    control = el_control(nthreads = 1)
  )
  fit2 <- el_glm(y ~ ., family = poisson("sqrt"), data = df)
  expect_equal(fit@optim, fit2@optim)
  expect_equal(fit@parTests, fit2@parTests)
  wfit <- el_glm(y ~ .,
    family = poisson("sqrt"), data = df, weights = w,
    control = el_control(nthreads = 1)
  )
  wfit2 <- el_glm(y ~ ., family = poisson("sqrt"), data = df, weights = w)
  expect_equal(wfit@optim, wfit2@optim)
  expect_equal(wfit@parTests, wfit2@parTests)
})

#' @srrstats {G5.4, G5.4a, G5.5} We perform correctness tests with two different
#'   data sets; one is the `women` data set from the `datasets` package and the
#'   other is simulated. `set.seed()` is used for the simulated data.
#'   Specifically, we consider the correctness of `elt()`, one of the main
#'   functions in the package. This function performs a constrained optimization
#'   of empirical likelihood. Since this is the first implementation of using
#'   projected gradient descent to the nested optimization of empirical
#'   likelihood, we consider a very simple problem with only two parameters.
#'   The hypothesis states that the second parameter in the model is some fixed
#'   number (150 for the `women` and 0 for the simulated data set). Then the
#'   optimization problem reduces to finding an optimal value for the first
#'   parameter, and this can be done using `optim()` function in the `stats`
#'   package. We test whether the `eld()` and the `optim()` produces the same
#'   results in terms of the solution and the convergence status. The tolerance
#'   for comparing floating point numbers is set to `1e-07`, which we used
#'   throughout all unit tests in the package. The threshold value `th` for
#'   empirical likelihood calculation is set to a large number `1e+10`. This
#'   prevents the optimization routine from exiting early when large empirical
#'   likelihood value is encountered. Note that it serves as a safeguard and
#'   does not affect the correctness of the test. No further manual
#'   specification of the control parameters are made in `el_control()` or
#'   `optim()`.
test_that("Correctness tests.", {
  skip_on_cran()
  fit <- el_mean(women, par = colMeans(women))
  out <- elt(fit, lhs = c(0, 1), rhs = 150, control = el_control(th = 1e+10))
  one_dim_function <- function(x) {
    -logLR(el_mean(women, par = c(x, 150)))
  }
  out2 <- optim(60, one_dim_function, method = "Brent", lower = 58, upper = 72)
  expect_true(out@optim$convergence)
  expect_identical(out2$convergence, 0L)
  expect_equal(unname(out@optim$par[1]), out2$par, tolerance = 1e-07)
  set.seed(65456)
  n <- 100
  x <- matrix(rnorm(n * 2), ncol = 2)
  fit2 <- el_mean(x, par = c(0, 0))
  out3 <- elt(fit2, lhs = c(0, 1), rhs = 0, control = el_control(th = 1e+10))
  one_dim_function2 <- function(y) {
    -logLR(el_mean(x, par = c(y, 0)))
  }
  out4 <- optim(-2, one_dim_function2, method = "Brent", lower = -4, upper = 4)
  expect_true(out3@optim$convergence)
  expect_identical(out4$convergence, 0L)
  expect_equal(out3@optim$par[1], out4$par, tolerance = 1e-07)
})

#' @srrstats {RE7.1, RE7.1a} We first simulate a design matrix with fixed
#'   coefficients. Then two response vectors `y1` and `y2` are generated, where
#'   `y1` has no noise while some noise is added to `y2`. Fitting a model using
#'   `el_lm()` is faster with `y1` than with `y2`.
test_that("Exact relationships between predictor and response.", {
  skip_on_cran()
  set.seed(311116)
  n <- 1000
  p <- 20
  b <- rep(1, p)
  x <- matrix(rnorm(n * p), ncol = p)
  y1 <- 1 + x %*% b
  y2 <- 1 + x %*% b + rnorm(n)
  df1 <- data.frame(y1, x)
  df2 <- data.frame(y2, x)
  out1 <- system.time(el_lm(y1 ~ ., df1))["elapsed"]
  out2 <- system.time(el_lm(y2 ~ ., df2))["elapsed"]
  expect_lte(out1, out2)
})

test_that("`el_glm()` (binomial - probit).", {
  skip_on_cran()
  set.seed(224543)
  n <- 100
  p <- 2
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 2), times = 50)
  l <- 1.5 + x %*% as.vector(b)
  mu <- pnorm(l)
  y <- vapply(mu, FUN = function(x) rbinom(1, 1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ ., family = binomial("probit"), data = df)
  wfit <- el_glm(y ~ ., family = binomial("probit"), data = df, weights = w)
  lhs <- list(
    matrix(c(1, -100, 0), nrow = 1),
    matrix(c(0, 1, -1), nrow = 1)
  )
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
  expect_s4_class(elmt(wfit, lhs = lhs), "ELMT")
})

test_that("`el_glm()` (binomial - log).", {
  skip_on_cran()
  set.seed(54324)
  n <- 100
  p <- 2
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 1.0005), times = 50)
  l <- -3 + x %*% as.vector(b)
  mu <- exp(l)
  y <- vapply(mu, FUN = function(x) rbinom(1, 1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ ., family = binomial("log"), data = df)
  wfit <- el_glm(y ~ ., family = binomial("log"), data = df, weights = w)
  lhs <- list(
    matrix(c(1, 5, 0), nrow = 1),
    matrix(c(0, 1, -1), nrow = 1)
  )
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
  expect_s4_class(elmt(wfit, lhs = lhs), "ELMT")
})

test_that("`el_glm()` (poisson - identity).", {
  skip_on_cran()
  set.seed(534)
  n <- 100
  p <- 3
  b <- rnorm(p, sd = 0.5)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 2), times = 50)
  l <- 3 + x %*% as.vector(b)
  mu <- l
  y <- vapply(mu, FUN = function(x) rpois(1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ ., family = poisson("identity"), data = df)
  wfit <- el_glm(y ~ ., family = poisson("identity"), data = df, weights = w)
  lhs <- list(
    matrix(c(1, 100, 0, 0), nrow = 1),
    matrix(c(0, 0, 1, -100), nrow = 1)
  )
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
  expect_s4_class(elmt(wfit, lhs = lhs), "ELMT")
})

test_that("`el_glm()` (poisson - sqrt).", {
  skip_on_cran()
  set.seed(15234)
  n <- 100
  p <- 3
  b <- rnorm(p, mean = 0, sd = .2)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- rep(c(1, 2), times = 50)
  l <- 0.5 + x %*% as.vector(b)
  mu <- l^2
  y <- vapply(mu, FUN = function(x) rpois(1, x), FUN.VALUE = integer(1))
  df <- data.frame(y, x)
  fit <- el_glm(y ~ ., family = poisson("sqrt"), data = df)
  wfit <- el_glm(y ~ ., family = poisson("sqrt"), data = df, weights = w)
  lhs <- list(
    matrix(c(1, -5, 0, 0), nrow = 1),
    matrix(c(0, 0, 1, -100), nrow = 1)
  )
  expect_s4_class(elmt(fit, lhs = lhs), "ELMT")
  expect_s4_class(elmt(wfit, lhs = lhs), "ELMT")
})

#' @srrstats {G5.9, G5.9b} Different random seeds do not produce significantly
#'   different critical values.
test_that("Noise susceptibility tests.", {
  skip_on_cran()
  fit <- el_lm(mpg ~ cyl + disp, data = mtcars)
  lhs <- list(matrix(c(1, 33, 0), nrow = 1), matrix(c(0, 1, -100), nrow = 1))
  set.seed(5246356)
  cv <- elmt(fit, lhs = lhs)@cv
  set.seed(195646)
  cv2 <- elmt(fit, lhs = lhs)@cv
  expect_equal(cv, cv2, tolerance = 1e-02)
})

test_that(
  "Parallel computation yields the same results (gaussian - inverse).",
  {
    skip_on_cran()
    fit <- el_glm(mpg ~ disp + hp + wt,
      family = gaussian("inverse"),
      data = mtcars, control = el_control(nthreads = 1)
    )
    fit2 <- el_glm(mpg ~ disp + hp + wt,
      family = gaussian("inverse"),
      data = mtcars
    )
    expect_equal(fit@optim, fit2@optim)
    expect_equal(fit@parTests, fit2@parTests)
    wfit <- el_glm(mpg ~ disp + hp + wt,
      family = gaussian("inverse"), data = mtcars, weights = gear,
      control = el_control(nthreads = 1)
    )
    wfit2 <- el_glm(mpg ~ disp + hp + wt,
      family = gaussian("inverse"),
      data = mtcars, weights = gear
    )
    expect_equal(wfit@optim, wfit2@optim)
    expect_equal(wfit@parTests, wfit2@parTests)
  }
)
