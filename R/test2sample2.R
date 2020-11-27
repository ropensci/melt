#' @export
## two sample test for mean
test2sample2 <- function(x, y, b = .1, maxit = 1000, abstol = 1e-8) {
  ### argument check(numeric, vector, not all same, etc...)

  ###
  result <- setNames(vector("list", 4),
                     c("par", "nlogLR", "iterations", "convergence"))

  # check convex hull constraint
  ub <- min(max(x), max(y))
  lb <- max(min(x), min(y))
  if (ub <= lb) {
    result$nlogLR <- Inf
    result$convergence <- -1
    return(result)
  }

  # initialization
  par <- (lb + ub) / 2
  nx <- length(x)
  ny <- length(y)
  N <- nx + ny
  alpha <- 1 / N
  convergence <- F
  iterations <- 0
  v <- 0

  # minimization
  while (convergence == F) {
    # lambda update
    #### separate function for lambda!!!
    lx <- el.mean(par, x)$lambda
    ly <- el.mean(par, y)$lambda

    # gradient
    grad <-
      c(sum(plog.prime(1 + lx * (x - par), threshold = 1 / nx)) * (-lx),
        sum(plog.prime(1 + ly * (y - par), threshold = 1 / ny)) * (-ly))

    # direction
    d <- dfp1d(grad)

    # direction change reverts momentum
    if (sign(d) != sign(v)) {
      v <- -v
    }

    # lb, ub update
    if (sign(d) > 0) {
      lb <- par
    } else {
      ub <- par
    }

    # convergence check & parameter update
    if (abs(d * sum(grad)) < abstol | ub - lb < abstol) {
      convergence = T
    } else {
      iterations <- iterations + 1
      # step halving to satisfy convex hull constraint
      v <- b * v + d
      while (par + alpha * v <= lb | par + alpha * v >= ub) {
        v <- v / 2
      }
      par <- par + alpha * v
      if (iterations == maxit)
        break
    }
  }

  # result
  result$par <- par
  result$nlogLR <- el.mean(par, x)$nlogLR + el.mean(par, y)$nlogLR
  result$iterations <- iterations
  result$convergence <- ifelse(convergence, 1, 0)
  result
}
