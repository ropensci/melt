#' @describeIn pVal Extracts the \eqn{p}-value.
setMethod("pVal", "EL", function(object, ...) {
  object@pval
})

#' @describeIn pVal Extracts the multiplicity adjusted \eqn{p}-values.
setMethod("pVal", "ELMT", function(object, ...) {
  object@pval
})

#' @describeIn pVal Extracts the \eqn{p}-value.
setMethod("pVal", "ELT", function(object, ...) {
  object@pval
})

#' @describeIn pVal Extracts the \eqn{p}-value.
setMethod("pVal", "SummaryEL", function(object, ...) {
  object@pval
})

#' @describeIn pVal Extracts the \eqn{p}-value.
setMethod("pVal", "SummaryELT", function(object, ...) {
  object@pval
})

#' @describeIn pVal Extracts the multiplicity adjusted \eqn{p}-values.
setMethod("pVal", "SummaryELMT", function(object, ...) {
  object@pval
})

#' @describeIn pVal Extracts the \eqn{p}-value.
setMethod("pVal", "SummaryLM", function(object, ...) {
  object@pval
})
