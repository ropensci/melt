#' @describeIn conv Extracts the convergence status of the model with respect to
#'   the Lagrange multiplier.
#' @srrstats {RE4.7} `conv()` method extracts information on convergence status.
setMethod("conv", "EL", function(object, ...) {
  getOptim(object)$convergence
})

#' @describeIn conv Extracts the convergence status of the model with respect to
#'   the parameter.
setMethod("conv", "CEL", function(object, ...) {
  getOptim(object)$convergence
})

#' @describeIn conv Extracts the convergence status of the model. See the
#'   documentation of \linkS4class{EL} and \linkS4class{CEL}.
setMethod("conv", "SummaryLM", function(object, ...) {
  object@convergence
})

#' @describeIn conv Extracts the convergence status of the model with respect to
#'   the parameter (or the Lagrange multiplier if `lhs` is `NULL`).
setMethod("conv", "ELT", function(object, ...) {
  getOptim(object)$convergence
})
