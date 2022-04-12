#' Empirical likelihood test for mean
#'
#' Computes empirical likelihood for mean parameter.
#'
#' @param par A numeric vector of parameter values to be tested.
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. Defaults to \code{NULL}, corresponding to identical
#'   weights. If non-\code{NULL}, weighted empirical likelihood is computed.
#' @param control A list of control parameters set by
#'   \code{\link{melt_control}}.
#' @param model A logical. If \code{TRUE} the model matrix used for fitting is
#'   returned.
#' @return A list of class \code{"el"} as described in \code{\link{el_eval}}.
#' @references Glenn, N.L., and Yichuan Zhao. 2007.
#'   “Weighted Empirical Likelihood Estimates and Their Robustness Properties.”
#'   Computational Statistics & Data Analysis 51 (10): 5130–41.
#'   \doi{10.1016/j.csda.2006.07.032}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1).
#'   \doi{10.1214/aos/1176347494}.
#' @seealso \link{el_eval}, \link{melt_control}
#' @examples
#' # scalar mean
#' par <- 0
#' x <- rnorm(100L)
#' el_mean(par, x)
#'
#' # vector mean
#' par <- c(0, 0)
#' x <- matrix(rnorm(100L), ncol = 2L)
#' el_mean(par, x)
#'
#' # weighted EL
#' par <- c(0, 0)
#' x <- matrix(rnorm(100), ncol = 2L)
#' w <- rep(c(1, 2), each = 25L)
#' el_mean(par, x, w)
#' @export
el_mean <- function(par, x, weights = NULL, control = melt_control(),
                    model = TRUE) {
  mm <- as.matrix(x)
  if (!is.numeric(mm) || !all(is.finite(mm)))
    stop("'x' must be a finite numeric matrix")
  if (nrow(mm) < 2L)
    stop("not enough 'x' observations")
  if (!is.numeric(par) || !all(is.finite(par)))
    stop("'par' must be a finite numeric vector")
  if (length(par) != ncol(mm))
    stop("'par' and 'x' have incompatible dimensions")
  if (!inherits(control, "melt_control") || !is.list(control))
    stop("invalid 'control' supplied")
  if (!is.null(weights)) {
    weights <- check_weights(weights, nrow(mm))
  }
  out <- mean_(par, mm, control$maxit_l, control$tol_l, control$th, weights)
  out$weights <- weights
  if (model)
    out$data.matrix <- mm
  class(out) <- "el"
  out
}
