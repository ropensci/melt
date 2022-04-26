#' @exportMethod eld
setGeneric("eld", function(object, control = control_el()) {
  standardGeneric("eld")
})

#' @exportMethod confreg
setGeneric("confreg", function(object, parm, level = 0.95,
                               cv = qchisq(level, 2L), npoints = 50L,
                               control = control_el()) {
  standardGeneric("confreg")
})

#' Plot methods
#'
#' Provides plot methods for objects that inherit from class \linkS4class{EL}.
#'
#' @param x Object to be plotted.
#' @param y Not used.
#' @param ... Further graphical parameters (see \code{\link[graphics]{par}}).
#' @usage NULL
#' @seealso \link{confreg}, \link{eld}
#' @exportMethod plot
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' Print methods
#'
#' Provides print methods for objects that inherit from class \linkS4class{EL}.
#'
#' @name print-method
#' @param x Object to be printed.
#' @param ... Further arguments passed to other methods.
#' @param digits Number of significant digits to be passed to
#'   \code{\link[base]{format}}.
#' @param signif.stars Logical. If \code{TRUE}, ‘significance stars’ are printed
#'   for each coefficient.
#' @usage NULL
#' @exportMethod print
setGeneric("print", function(x, ...) standardGeneric("print"))

#' Summary methods
#'
#' Provides summary methods for objects that inherit from class
#'   \linkS4class{EL}.
#'
#' @name summary-method
#' @param object Object for which a summary is desired.
#' @param ... Additional arguments affecting the summary produced.
#' @usage NULL
#' @exportMethod summary
setGeneric("summary", function(object, ...) standardGeneric("summary"))

