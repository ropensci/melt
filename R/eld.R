#' Empirical likelihood displacement
#'
#' Empirical likelihood displacement for objects inheriting from class
#'   \code{"el_test"}.
#'
#' @param object A fitted \code{"el_test"} object.
#' @param control A list of control parameters. See ‘Details’.
#' @export
eld <- function(object, control = list()) {
  optcfg <- check_control(control)
  maxit <- optcfg$maxit
  tol <- optcfg$tol
  th <- optcfg$th
  w <- object$weights
  out <- eld_(object$optim$method, object$coefficients, object$data.matrix,
              maxit, tol, th , w)
  out
}
