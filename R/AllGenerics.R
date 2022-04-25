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


#' Plot methods for objects that inherit from class \linkS4class{EL}
#'
#' Provides plot methods for objects of class \linkS4class{ConfregEL} or
#'   \linkS4class{ELD}.
#'
#' @param x Object for plotting.
#' @param y Not used.
#' @param ... Further graphical parameters (see \code{\link[graphics]{par}}).
#' @usage NULL
#' @seealso \link{confreg}, \link{eld}
#' @exportMethod plot
setGeneric("plot", function(x, y, ...) standardGeneric("plot"), signature = "x")


#' #' Plot methods for objects that inherit from class \linkS4class{EL}
#' #'
#' #' Provides plot methods for S4 objects that inherit from class \linkS4class{EL}.
#' #'
#' #' @param x An S4 object of class \code{\link{EL}}.
#' #' @param y An S4 object of class \code{\link{EL}}.
#' #' @param digits An S4 object of class \code{\link{EL}}.
#' #' @param ... further arguments passed to or from other methods.
#' setGeneric("print", function(x, y, digits, ...) standardGeneric("print"))
#'
#'
#' setMethod("show", "EL", function(object) print(object))
