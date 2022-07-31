#' @rdname getDF
setMethod("getDF", "EL", function(object) {
  object@df
})

#' @rdname getDF
setMethod("getDF", "logLikEL", function(object) {
  object@df
})

#' @rdname getDF
setMethod("getDF", "SummaryLM", function(object) {
  object@df
})
