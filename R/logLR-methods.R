#' @rdname logLR
setMethod("logLR", "EL", function(object, ...) {
  object@loglr
})

#' @rdname logLR
setMethod("logLR", "ELT", function(object, ...) {
  object@loglr
})

#' @rdname logLR
setMethod("logLR", "SummaryEL", function(object, ...) {
  object@loglr
})

#' @rdname logLR
setMethod("logLR", "SummaryELT", function(object, ...) {
  object@loglr
})

#' @rdname logLR
setMethod("logLR", "SummaryLM", function(object, ...) {
  object@loglr
})
