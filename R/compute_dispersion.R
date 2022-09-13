#' Compute dispersion
#'
#' Compute dispersion for an object of class \linkS4class{QGLM} in [elt()].
#'
#' @param object An object of class \linkS4class{QGLM}.
#' @param par A numeric vector of the given parameter value.
#' @details [compute_dispersion()] is called only when the `"boot"` option in
#'   the `calibrate` argument is selected. It adjusts the dispersion parameter
#'   so that the null transformation can be applied to the data without
#'   contradiction.
#' @return A single numeric for the dispersion with the given parameter value.
#' @noRd
compute_dispersion <- function(object, par) {
  mm <- getData(object)
  n <- nobs(object)
  s <- mm[, 1L]
  y <- mm[, 2L]
  x <- mm[, -c(1L, 2L)]
  if (is.null(dim(x))) {
    dim(x) <- c(n, ncol(mm) - 2L)
  }
  yhat <- object@family$linkinv(x %*% par + s)
  v <- object@family$variance(yhat)
  w <- weights(object)
  if (is.null(w)) {
    sum((y - yhat)^2L / v) / n
  } else {
    sum(w * ((y - yhat)^2L / v)) / n
  }
}
