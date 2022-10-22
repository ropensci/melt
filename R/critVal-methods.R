#' @rdname critVal
setMethod("critVal", "ELMT", function(object, ...) {
  object@cv
})

#' @rdname critVal
setMethod("critVal", "ELT", function(object, ...) {
  object@cv
})

#' @rdname critVal
setMethod("critVal", "SummaryELMT", function(object, ...) {
  object@cv
})

#' @rdname critVal
setMethod("critVal", "SummaryELT", function(object, ...) {
  object@cv
})
