#' @rdname logProb
setMethod("logProb", "EL", function(object, ...) {
  object@logp
})

#' @rdname logProb
setMethod("logProb", "ELT", function(object, ...) {
  object@logp
})
