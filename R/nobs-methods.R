#' @describeIn EL Extracts the number of observations from a model.
#' @exportMethod nobs
#' @usage NULL
#' @srrstats {RE4.5} `nobs()` method extracts the number of observations in the
#'   model.
setMethod("nobs", "EL", function(object, ...) {
  object@nobs
})
