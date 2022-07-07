#' @describeIn conv Extracts the convergence status of the model with respect to
#'   the Lagrange multiplier.
setMethod("conv", "EL", function(object, ...) {
  object@optim$convergence
})

#' @describeIn conv Extracts the convergence status of the model with respect to
#'   the parameter.
setMethod("conv", "CEL", function(object, ...) {
  object@optim$convergence
})

#' @describeIn conv Extracts the convergence status of the model with respect to
#'   the parameter (or the Lagrange multiplier if `lhs` is `NULL`).
setMethod("conv", "ELT", function(object, ...) {
  object@optim$convergence
})
