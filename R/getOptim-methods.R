#' @rdname getOptim
setMethod("getOptim", "EL", function(object, ...) {
  object@optim
})

#' @rdname getOptim
setMethod("getOptim", "ELT", function(object, ...) {
  object@optim
})
