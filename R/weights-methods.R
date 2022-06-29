#' @rdname weights
setMethod("weights", "EL", function(object, ...) {
  if (!missing(...)) {
    warning("extra arguments are not supported.")
  }
  out <- getWeights(object)
  if (length(out) == 0L) {
    return(NULL)
  }
  out
})
