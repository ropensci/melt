#' Empirical likelihood test for mean
#'
#' Computes empirical likelihood for mean parameter.
#'
#' @param par Numeric vector of parameters to be tested.
#' @param x Numeric matrix or vector of data. For matrix \code{x}, each row corresponds to an observation.
#' @param control A list of control parameters. See ‘Details’.
#' @return A list with class \code{"el_test"}.
#'
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.” The Annals of Statistics 18 (1). \doi{10.1214/aos/1176347494}.
#' @examples
#' ## scalar mean
#' par <- 0
#' x <- rnorm(100)
#' el_mean(par, x)
#'
#' ## vector mean
#' x <- matrix(rnorm(100), ncol = 2)
#' par <- c(0, 0)
#' el_mean(par, x)
#' @export
el_mean <- function(par, x, control = list()) {
  if (!is.numeric(par)) stop("'par' must be a numeric vector")
  optcfg <- check_control(control)
  out <- EL_mean(par, x, optcfg$maxit, optcfg$abstol, optcfg$threshold)
  out$data <- x
  out$data.name <- deparse1(substitute(x))
  return(out)
}
