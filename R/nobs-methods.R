#' @rdname nobs
setMethod("nobs", "EL", function(object, ...) {
  object@nobs
})

#' @rdname nobs
setMethod("nobs", "SummaryEL", function(object, ...) {
  object@nobs
})

#' @rdname nobs
setMethod("nobs", "SummaryLM", function(object, ...) {
  object@nobs
})
