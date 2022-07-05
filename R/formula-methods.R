#' @describeIn LM Extracts the symbolic model formula used in [`el_lm()`] or
#'   [`el_glm()`].
#' @exportMethod formula
#' @usage NULL
#' @srrstats {RE4.4} `formula()` method extracts the model formula used.
setMethod("formula", "LM", function(x, ...) {
  formula(x@terms)
})
