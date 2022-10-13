#' @rdname logL
setMethod("logL", "EL", function(object, ...) {
  object@logl
})

#' @rdname logL
setMethod("logL", "ELT", function(object, ...) {
  object@logl
})

#' @rdname logL
setMethod("logL", "SummaryEL", function(object, ...) {
  object@logl
})

#' @rdname logL
setMethod("logL", "SummaryELT", function(object, ...) {
  object@logl
})

#' @rdname logL
setMethod("logL", "SummaryLM", function(object, ...) {
  object@logl
})
