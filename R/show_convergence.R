#' Shows a convergence status
#'
#' Shows a convergence status for model objects.
#'
#' @param x An object to be printed.
#' @return A single character vector for the convergence status.
#' @noRd
show_convergence <- function(x) {
  iter <- getOptim(x)$iterations
  cstr <- getOptim(x)$cstr
  if (cstr) {
    type <- "Constrained EL:"
    status <- if (conv(x)) {
      "converged"
    } else {
      if (iter == getControlEL(x)@maxit) {
        "maximum iterations reached"
      } else if (iter == 0L) {
        "initialization failed"
      } else {
        "not converged"
      }
    }
  } else {
    type <- "EL evaluation:"
    status <- if (conv(x)) {
      "converged"
    } else {
      if (iter == getControlEL(x)@maxit_l) {
        "maximum iterations reached"
      } else {
        "solver error"
      }
    }
  }
  paste(type, status)
}
