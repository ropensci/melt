#' An S4 class \code{"ConfregEL"}
#'
#' An S4 class for a two-dimensional confidence region.
#'
#' @aliases  ConfregEL
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
