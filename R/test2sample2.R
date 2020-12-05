#' Two sample test for equal mean
#'
#' Two sample test for equal mean
#'
#' @param x a vector of data for one  group.
#' @param y a vector of data for the other  group.
#' @param b a momentum parameter for minimization. Defaults to .9.
#' @param alpha an optional step size. Defaults to 1.
#' @param maxit an optional value for the maximum number of iterations. Defaults to 1000.
#' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' test2sample2(x, y)
#'
#' @export
test2sample2 <- function(x, y, b = 0.9, alpha = 1,
                         maxit = 1000,  abstol = 1e-8) {
  result <- vector("list", 4)
  names(result) <- c("par", "nlogLR", "iterations", "convergence")

  # check convex hull constraint
  xbar <- mean(x)
  ybar <- mean(y)
  ub <- min(max(x), max(y), max(xbar, ybar))
  lb <- max(min(x), min(y), min(xbar, ybar))
  if (ub <= lb) {
    result$nlogLR <- Inf
    result$convergence <- -1
    return(result)
  }

  # initialization
  par <- (lb + ub) / 2
  nx <- length(x)
  ny <- length(y)
  N <- max(nx, ny)
  convergence <- 0
  iterations <- 0
  v <- 0

  # minimization
  while (convergence == 0) {
    # lambda update
    lx <- el.mean(par, x)$lambda
    ly <- el.mean(par, y)$lambda
    # gradient
    grad <-
      c(sum(plog.prime(1 + lx * (x - par), threshold = 1 / nx)) * (-lx),
        sum(plog.prime(1 + ly * (y - par), threshold = 1 / ny)) * (-ly)) / N
    # direction
    d <- dfp1d(grad)
    # direction change reverts momentum
    if (sign(d) != sign(v)) {
      v <- 0
      alpha <- alpha / 2
    }
    # lb, ub update
    if (sign(d) > 0) {
      lb <- par
    } else {
      ub <- par
    }
    # convergence check & parameter update
    if ((abs(d * sum(grad)) < abstol | ub - lb < abstol) & iterations > 0) {
      convergence <- 1
    } else {
      iterations <- iterations + 1
      # step halving to satisfy convex hull constraint
      v <- b * v + d
      while (par + alpha * v <= lb | par + alpha * v >= ub) {
        alpha <- alpha / 2
      }
      par <- par + alpha * v
      if (iterations == maxit)
        break
    }
  }

  # result
  result$par <- par
  result$nlogLR <- sum(log(1 + (x - par) * lx)) + sum(log(1 + (y - par) * ly))
  result$iterations <- iterations
  result$convergence <- convergence
  result
}
