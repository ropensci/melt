#' Empirical likelihood test for mean
#'
#' Computes empirical likelihood for mean parameter.
#'
#' @param par A numeric vector of parameter values to be tested.
#' @param x A numeric matrix, or an object that can be coerced to a numeric matrix. Each row corresponds to an observation.
#' @param weights An optional numeric vector of weights. Defaults to \code{NULL}, corresponding to identical weights. If non-\code{NULL}, weighted empirical likiehood is computed.
#' @param control A list of control parameters. See ‘Details’.
#' @return A list with class \code{"el_test"}.
#'
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.” The Annals of Statistics 18 (1). \doi{10.1214/aos/1176347494}.
#' @references Glenn, N.L., and Yichuan Zhao. 2007. “Weighted Empirical Likelihood Estimates and Their Robustness Properties.” Computational Statistics & Data Analysis 51 (10): 5130–41. \doi{10.1016/j.csda.2006.07.032}.
#' @examples
#' ## scalar mean
#' par <- 0
#' x <- rnorm(100)
#' el_mean(par, x)
#'
#' ## vector mean
#' par <- c(0, 0)
#' x <- matrix(rnorm(100L), ncol = 2)
#' el_mean(par, x)

#' ## weighted EL
#' par <- c(0, 0)
#' x <- matrix(rnorm(100), ncol = 2)
#' el_mean(par, x, weights = rep(c(1,2), each = 25))
#' @export
el_mean <- function(par, x, weights = NULL, control = list())
{
  if (!is.numeric(par) || any(!is.finite(par)))
    stop("'par' must be a finite numeric vector")
  mm <- as.matrix(x)
  if (!is.numeric(mm) || any(!is.finite(mm)))
    stop("'x' must be a finite numeric matrix")
  if (NROW(mm) < 2L)
    stop("not enough 'x' observations")
  if (length(par) != NCOL(mm))
    stop("'par' and 'x' have incompatible dimensions")
  if (is.null(weights)) {
    w <- rep(1, NROW(mm))
  } else {
    w <- as.numeric(weights)
  }
  if (any(!is.finite(w)))
    stop("'weights' must be a finite numeric vector")
  if (any(w < 0))
    stop("negative 'weights' are not allowed")
  if (length(w) != NROW(mm))
    stop("'x' and 'weights' have incompatible dimensions")
  optcfg <- check_control(control)
  # out <- EL_mean(par, mm, optcfg$maxit, optcfg$abstol, optcfg$threshold)
  out <- EL_mean_weight(par, mm, w,
                        optcfg$maxit, optcfg$abstol, optcfg$threshold)
  out$data.matrix <- mm
  out$data.name <- deparse1(substitute(x))
  return(out)
}
