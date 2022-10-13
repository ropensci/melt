#' @rdname getOptim
setMethod("getOptim", "EL", function(object, ...) {
  object@optim
})

#' @rdname getOptim
setMethod("getOptim", "ELT", function(object, ...) {
  object@optim
})

#' @rdname getOptim
setMethod("getOptim", "SummaryEL", function(object, ...) {
  object@optim
})

#' @rdname getOptim
setMethod("getOptim", "SummaryELT", function(object, ...) {
  object@optim
})

#' @rdname getOptim
setMethod("getOptim", "SummaryLM", function(object, ...) {
  object@optim
})
