#' Empirical log-likelihood
#'
#' Extracts empirical log-likelihood of a model represented evaluated at the
#'   estimated coefficients. Package \strong{melt} adds a method for objects
#'   inheriting from class \linkS4class{EL}.
#'
#' @param object A fitted \linkS4class{EL} object.
#' @param ... Some methods for this generic function require extra arguments.
#'   None are used in this method.
#' @return An object of class \code{"logLik"} with an attribute \code{df} that
#'   gives the number of (estimated) parameters in the model.
#' @examples
#' fit <- el_lm(formula = mpg ~ wt, data = mtcars)
#' logLik(fit)
setMethod(
  "logLik", "EL",
  function(object, ...) {
    if (!missing(...)) {
      warning("extra arguments are not supported")
    }
    p <- object@npar
    rhs <- object@coefficients
    out <- lht(object, rhs = rhs)
    val <- out@logl
    attr(val, "df") <- p
    val
  }
)
