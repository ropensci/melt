#' @describeIn coef Extracts the numeric vector of the maximum empirical
#'   likelihood estimates.
setMethod("coef", "EL", function(object, ...) {
  object@coefficients
})

#' @describeIn coef Extracts the list of numeric vectors of the maximum
#'   empirical likelihood estimates. Each element of the list corresponds to a
#'   distinct hypothesis.
setMethod("coef", "ELMT", function(object, ...) {
  object@coefficients
})
