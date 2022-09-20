#' @describeIn getDF Extracts the degrees of freedom.
setMethod("getDF", "EL", function(object) {
  object@df
})

#' @describeIn getDF Extracts the vector of marginal degrees of freedom.
setMethod("getDF", "ELMT", function(object) {
  object@df
})

#' @describeIn getDF Extracts the (chi-square) degrees of freedom.
setMethod("getDF", "ELT", function(object) {
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
