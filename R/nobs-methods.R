#' @rdname nobs
#' @srrstats {RE4.5} `nobs()` method extracts the number of observations from a
#'   model.
setMethod("nobs", "EL", function(object, ...) {
  object@nobs
})
