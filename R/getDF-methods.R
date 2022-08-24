#' @describeIn getDF Extracts the degrees of freedom.
setMethod("getDF", "EL", function(object) {
  object@df
})

#' @describeIn getDF Extracts the marginal degrees of freedoms.
setMethod("getDF", "ELMT", function(object) {
  object@df
})

#' @describeIn getDF Extracts the degrees of freedom.
setMethod("getDF", "logLikEL", function(object) {
  object@df
})

#' @describeIn getDF Extracts the degrees of freedom.
setMethod("getDF", "SummaryLM", function(object) {
  object@df
})
