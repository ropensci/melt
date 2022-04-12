check_weights <- function(weights, nw) {
  if (!is.numeric(weights))
    stop("'weights' must be a numeric vector")
  w <- as.numeric(weights)
  if (!all(is.finite(w)))
    stop("'weights' must be a finite numeric vector")
  if (any(w < 0))
    stop("negative 'weights' not allowed")
  if (length(w) != nw)
    stop("length of 'weights' is incompatible with data")
  w <- (nw / sum(w)) * w
  w
}

check_rhs <- function(rhs, p) {
  rhs <- as.vector(rhs, "numeric")

  if (!is.numeric(rhs) || !all(is.finite(rhs)))
    stop("'rhs' must be a finite numeric vector")
  if (length(rhs) != p)
    stop(paste("length of 'rhs' should be "), p)
  rhs
}

check_lhs <- function(lhs, p) {
  lhs <- as.matrix(lhs)

  if (!is.numeric(lhs) || !all(is.finite(lhs)))
    stop("'lhs' must be a finite numeric matrix")
  if (ncol(lhs) != p)
    stop("'object' and 'lhs' have incompatible dimensions")
  q <- nrow(lhs)
  if (q == 0L || q > p)
    stop("'object' and 'lhs' have incompatible dimensions")
  lhs
}

check_hypothesis <- function(lhs, rhs, p) {
  if (is.null(rhs) && is.null(lhs)) {
    stop("either 'rhs' or 'lhs' must be provided")
  } else if (is.null(lhs)) {
    rhs <- check_rhs(rhs, p)
  } else if (is.null(rhs)) {
    lhs <- check_lhs(lhs, p)
    rhs <- rep(0, nrow(lhs))
  } else {
    lhs <- check_lhs(lhs, p)
    rhs <- check_rhs(rhs, nrow(lhs))
  }
  list(l = lhs, r = rhs)
}
