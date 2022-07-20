#' @rdname logL
setMethod("logL", "EL", function(object, ...) {
  object@logl
})

#' @rdname logL
setMethod("logL", "ELT", function(object, ...) {
  object@logl
})
