#' @rdname logLR
setMethod("logLR", "EL", function(object, ...) {
  object@loglr
})

#' @rdname logLR
setMethod("logLR", "ELT", function(object, ...) {
  object@loglr
})
