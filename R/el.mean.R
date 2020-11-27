#' Compute empirical likelihood for mean
#'
#' Compute empirical likelihood for mean
#'
#' @param theta a vector of parameters to be tested.
#' @param x a matrix or vector of data. Each row is an observation vector.
#' @param abstol an optional value for the absolute convergence tolerance. Defaults to 1e-8.
#' @param maxit an optional value for the maximum number of iterations. Defaults to 50.
#'
#' @return An object with S3 class "el"
#'
#' @examples
#' theta <- 0
#' x <- rnorm(100)
#' el.mean(theta, x)
#'
#' @import Matrix
#' @export
el.mean <- function(theta, x, abstol = 1e-8, maxit = 50) {
  ## Data
  x <- as.matrix(x)
  if (length(theta) != ncol(x))
    stop("length(theta) != ncol(x)")
  if (rankMatrix(x) != ncol(x))
    stop("model matrix must have full column rank")
  colnames(x) <- NULL
  n <- nrow(x)
  p <- ncol(x)
  # centering(estimating equation with mean-zero)
  z <- x - matrix(rep(theta, n), nrow = n, byrow = T)


  ## Minimization
  lambda <- rep(0, p)
  for (i in 1:maxit) {
    # function evaluation(initial)
    f0 <- -sum(plog(1 + z %*% lambda, threshold = 1 / n))
    # J matrix for least square
    J <- sqrt(-plog.dprime(1 + z %*% lambda, threshold = 1 / n)) * z
    # Y vector for least square
    Y <- plog.prime(1 + z %*% lambda, threshold = 1 / n) /
      sqrt(-plog.dprime(1 + z %*% lambda, threshold = 1 / n))
    # update lambda by NR method with least square & step halving
    r <- 1
    lambda_c <- lambda + solve(crossprod(J), crossprod(J, Y)) / r
    f1_c <- -sum(plog(1 + z %*% lambda_c, threshold = 1 / n))
    while (f1_c > f0) {
      r <- 2 * r
      lambda_c <- lambda + solve(crossprod(J), crossprod(J, Y)) / r
      f1_c <- -sum(plog(1 + z %*% lambda_c, threshold = 1 / n))
    }
    lambda <- lambda_c
    # function evaluation(update)
    f1 <- f1_c
    # convergence check
    if (f0 - f1 < abstol) {
      iterations <- i
      break
    }
    # increase iteration number
    if (i == maxit) {
      iterations <- maxit
      break
    } else {
      i <- i + 1
    }
  }


  ## Result
  result <- list()
  class(result) <- "el"
  result$nlogLR <- sum(log(1 + z %*% lambda))
  result$w <- drop((1 + z %*% lambda)^(-1) / n)
  result$lambda <- drop(lambda)
  result$grad <- -colSums(plog.prime(1 + z %*% lambda, threshold = 1 / n) * z)
  result$hessian <- -crossprod(z, plog.dprime(1 + z %*% lambda, threshold = 1 / n) * z)
  result$iterations <- iterations
  result$message <- ifelse(all.equal(sum(result$w), 1) == T,
                           "convex hull constraint satisfied",
                           "something wrong"
  )
  result
}
