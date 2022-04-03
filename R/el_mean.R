#' Empirical likelihood test for mean
#'
#' Computes empirical likelihood for mean parameter.
#'
#' @param par A numeric vector of parameter values to be tested.
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. If not provided, identical weights are applied. Otherwise,
#'   weighted empirical likelihood is computed.
#' @param control A list of control parameters. See ‘Details’ in
#'   \code{\link{el_eval}}.
#' @param model A logical. If \code{TRUE} the model matrix used for fitting is
#'   returned.
#' @return A list with class \code{"el_test"} as described in
#'   \code{\link{lht}}.
#' @references Glenn, N.L., and Yichuan Zhao. 2007.
#'   “Weighted Empirical Likelihood Estimates and Their Robustness Properties.”
#'   Computational Statistics & Data Analysis 51 (10): 5130–41.
#'   \doi{10.1016/j.csda.2006.07.032}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1).
#'   \doi{10.1214/aos/1176347494}.
#' @seealso \link{el_eval}, \link{lht}
#' @examples
#' # scalar mean
#' par <- 0
#' x <- rnorm(100)
#' el_mean(par, x)
#'
#' # vector mean
#' par <- c(0, 0)
#' x <- matrix(rnorm(100L), ncol = 2)
#' el_mean(par, x)
#'
#' # weighted EL
#' par <- c(0, 0)
#' x <- matrix(rnorm(100), ncol = 2)
#' w <- rep(c(1, 2), each = 25)
#' el_mean(par, x, w)
#' @export
el_mean <- function(par, x, weights, control = list(), model = TRUE) {
  mm <- as.matrix(x)
  if (!is.numeric(mm) || !all(is.finite(mm)))
    stop("'x' must be a finite numeric matrix")
  if (nrow(mm) < 2L)
    stop("not enough 'x' observations")
  if (!is.numeric(par) || !all(is.finite(par)))
    stop("'par' must be a finite numeric vector")
  if (length(par) != ncol(mm))
    stop("'par' and 'x' have incompatible dimensions")

  # check control
  optcfg <- check_control(control)
  if (missing(weights)) {
    out <- mean_(par, mm, optcfg$maxit, optcfg$tol, optcfg$th)
  } else {
    w <- check_weights(weights, nrow(mm))
    out <- mean_w_(par, mm, w, optcfg$maxit, optcfg$tol, optcfg$th)
    out$weights <- w
  }
  if (model)
    out$data.matrix <- mm
  class(out) <- "el_test"
  out
}
