#' @export
el_eval <- function(g, weights = NULL, control = list()) {
  # check g

  # check control
  optcfg <- check_control(control)
  if (is.null(weights)) {
    out <- EL_eval(g, optcfg$maxit, optcfg$abstol, optcfg$threshold)
  } else {
    if (!is.numeric(weights))
      stop("'weights' must be a numeric vector")
    w <- as.numeric(weights)
    if (any(!is.finite(w)))
      stop("'weights' must be a finite numeric vector")
    if (any(w < 0))
      stop("negative 'weights' are not allowed")
    if (length(w) != NROW(g))
      stop("'g' and 'weights' have incompatible dimensions")
    w <- (NROW(g) / sum(w)) * w
    out <- WEL_eval(g, w, optcfg$maxit, optcfg$abstol, optcfg$threshold)
  }
  out
}

g_lm <- function(y, x, par) {
  x * as.vector(y - x %*% par)
}
