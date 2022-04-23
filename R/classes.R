#' An S4 class \code{ConfregEL}
#'
#' An S4 class for a two-dimensional confidence region.
#'
#' @aliases ConfregEL
#' @slot points A numeric matrix with two columns for boundary points of a
#'   confidence region.
#' @slot estimates A numeric vector of length two for parameter estimates.
#' @slot level A confidence level required.
#' @slot cv A critical value for calibration of empirical likelihood ratio
#'   statistic.
#' @slot pnames A character vector of length two for the name of parameters.
setClass("ConfregEL",
  slots = c(
    points = "matrix", estimates = "numeric", level = "numeric", cv = "numeric",
    pnames = "character"
  ),
  prototype = list(
    points = NULL, estimates = NA_real_, level = NA_real_, cv = NA_real_,
    pnames = NA_character_
  )
)

#' An S4 class \code{ELD}
#'
#' An S4 class for empirical likelihood displacement.
#'
#' @aliases ELD
#' @slot eld A numeric vector of empirical likelihood displacement values.
setClass("ELD",
  slots = c(eld = "numeric"),
  prototype = list(eld = NA_real_)
)

#' Plot methods for package \strong{melt}
#'
#' Provides plot methods for S4 objects that inherit from class \code{"el"}.
#'
#' @param x An S4 object of class \code{\link{ConfregEL}} or \code{\link{ELD}}.
#' @param y y From the generic \code{plot} function, ignored for class objects.
#' @param ... A character vector of length two for the name of parameters.
#' @name plot-method
#' @rdname plot-method
#' @aliases plot
#' @exportMethod plot
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' @seealso \link{confreg}, \link{eld}
#' @examples
#' data("clothianidin")
#' fit <- el_lm(clo ~ trt, clothianidin)
#' eld <- eld(fit)
#' plot(eld)
#' @rdname plot-method
#' @importFrom graphics polygon
setMethod("plot", "ConfregEL", function(x, y, tt = NULL, ...) {
  plot(x@points[, 1L], x@points[, 2L], ...)
  polygon(x@points[, 1L], x@points[, 2L], ...)
})
#' @rdname plot-method
setMethod("plot", "ELD", function(x, y, ...,
                                  main = "Empirical Likelihood Displacement",
                                  ylab = "ELD", pch = 21) {
  plot(x@eld, ..., main = main, ylab = ylab, pch = pch)
})



#' #' Plot for eld objects
#' #'
#' #' Takes a fitted \code{"eld"} object and plots the empirical likelihood
#' #'   displacement values versus the observation index.
#' #'
#' #' @param x An object of class \code{"eld"}.
#' #' @param ... Additional arguments to be passed to \code{\link[base]{plot}}.
#' #' @param main A title for the plot.
#' #' @param ylab A label for the y-axis.
#' #' @param pch A vector of plotting characters or symbols to use.
