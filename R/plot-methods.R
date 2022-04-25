#' @describeIn plot Plots a two-dimensional confidence region for model
#'   parameters.
#' @importFrom graphics polygon
setMethod("plot", "ConfregEL", function(x, y, ...) {
  plot(x@points[, 1L], x@points[, 2L], ...)
  polygon(x@points[, 1L], x@points[, 2L], ...)
})

#' @describeIn plot Plots empirical likelihood displacement values versus
#'   observation index.
#' @exportMethod plot
setMethod("plot", "ELD", function(x, y, ...) {
  args <- list(...)
  if (!exists("xlab", args)) {
    args$main <- "Empirical Likelihood Displacement"
  }
  if (!exists("ylab", args)) {
    args$ylab <- "ELD"
  }
  if (!exists("pch", args)) {
    args$pch <- 21
  }
  args$x <- x@eld
  do.call(plot, args)
})
