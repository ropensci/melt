#' @describeIn plot Plots a two-dimensional confidence region for model
#'   parameters.
setMethod("plot", "ConfregEL", function(x, ...) {
  args <- list(...)
  if (!exists("main", args)) {
    args$main <- paste(
      if (!is.na(x@level)) sprintf("%0.1f%%", 100 * x@level),
      "Confidence Region"
    )
  }
  if (!exists("xlab", args)) {
    args$xlab <- x@pnames[1L]
  }
  if (!exists("ylab", args)) {
    args$ylab <- x@pnames[2L]
  }
  if (!exists("pch", args)) {
    args$pch <- 20L
  }
  if (!exists("cex", args)) {
    args$cex <- 0.5
  }
  args$x <- getDataPart(x)[, 1L]
  args$y <- getDataPart(x)[, 2L]
  do.call(plot.default, args)
  do.call(polygon, args)
  points(x@estimates[1L], x@estimates[2L], pch = args$pch, cex = args$cex)
  text(x@estimates[1L], x@estimates[2L], expression(hat(theta)), pos = 4L)
})

#' @describeIn plot Plots empirical likelihood displacement values versus
#'   observation index. `eld()` is called implicitly.
setMethod("plot", "EL", function(x, ...) {
  plot(eld(x), ...)
})

#' @describeIn plot Plots empirical likelihood displacement values versus
#'   observation index.
setMethod("plot", "ELD", function(x, ...) {
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
  args$x <- getDataPart(x)
  do.call(plot.default, args)
})
