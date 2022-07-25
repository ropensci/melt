#' @rdname getSigTests
setMethod("getSigTests", "LM", function(object) {
  object@sigTests
})
