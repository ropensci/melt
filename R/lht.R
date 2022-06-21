#' Empirical likelihood test
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' \code{\link{lht}} was renamed to \code{\link{elt}} to avoid conflicts with
#' other packages.
#' @seealso \link{elt}
#' @keywords internal
#' @export
lht <- function(object, rhs = NULL, lhs = NULL, alpha = 0.05,
                calibrate = c("chisq", "boot", "f"), control = el_control()) {
  lifecycle::deprecate_warn("1.5.3", "lht()", "elt()")
  elt(object, rhs, lhs, alpha, calibrate, control)
}
