#' Empirical likelihood for the mean
#'
#' Computes empirical likelihood for the mean.
#'
#' @param par A numeric vector of parameter values to be tested.
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. Defaults to \code{NULL}, corresponding to identical
#'   weights. If non-\code{NULL}, weighted empirical likelihood is computed.
#' @param control A list of control parameters set by
#'   \code{\link{control_el}}.
#' @param model A logical. If \code{TRUE} the data matrix used for fitting is
#'   returned.
#' @return A list of class \code{"el"} with the following components:
#'   \item{optim}{A list with the following optimization results:
#'     \itemize{
#'       \item{\code{method } }{A character for method dispatch in internal
#'       functions.}
#'       \item{\code{par } }{The value at which empirical likelihood is
#'       evaluated.}
#'       \item{\code{lambda } }{The Lagrange multiplier of dual problem.}
#'       \item{\code{logLR } }{The (weighted) empirical log-likelihood ratio
#'       value.}
#'       \item{\code{iterations } }{The number of iterations performed.}
#'       \item{\code{convergence } }{A logical vector. \code{TRUE} indicates
#'       convergence of the algorithm.}
#'     }
#'   }
#'   \item{log.prob}{The log probabilities.}
#'   \item{loglik}{The log likelihood value.}
#'   \item{statistic}{The chi-square statistic.}
#'   \item{df}{The degrees of freedom of the statistic.}
#'   \item{p.value}{The \eqn{p}-value of the statistic.}
#'   \item{npar}{The number of parameters.}
#'   \item{weights}{The rescaled weights if non-\code{NULL} \code{weights} is
#'   supplied}
#'   \item{data.matrix}{The data matrix used for fitting if \code{model} is
#'   \code{TRUE}.}
#'   \item{coefficients}{The maximum empirical likelihood estimates of the
#'   parameters.}
#' @references Glenn, N.L., and Yichuan Zhao. 2007.
#'   “Weighted Empirical Likelihood Estimates and Their Robustness Properties.”
#'   Computational Statistics & Data Analysis 51 (10): 5130–41.
#'   \doi{10.1016/j.csda.2006.07.032}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{control_el}, \link{el_eval}, \link{lht}
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
#' # weighted data
#' par <- c(0, 0)
#' x <- matrix(rnorm(100), ncol = 2L)
#' w <- rep(c(1, 2), each = 25L)
#' el_mean(par, x, w)
#' @export
el_mean <- function(par, x, weights = NULL, control = control_el(),
                    model = TRUE) {
  mm <- as.matrix(x)
  n <- nrow(mm)
  p <- ncol(mm)
  if (!is.numeric(mm) || !all(is.finite(mm))) {
    stop("'x' must be a finite numeric matrix")
  }
  if (n < 2L) {
    stop("not enough 'x' observations")
  }
  if (get_rank_(mm) != p || n <= p) {
    stop("'x' must have full column rank")
  }
  if (!is.numeric(par) || !all(is.finite(par))) {
    stop("'par' must be a finite numeric vector")
  }
  if (length(par) != p) {
    stop("'par' and 'x' have incompatible dimensions")
  }
  if (!inherits(control, "control_el") || !is.list(control)) {
    stop("invalid 'control' supplied")
  }
  w <- check_weights(weights, n)
  if (!is.null(weights)) {
    cf <- colSums(mm * w) / n
  } else {
    cf <- colMeans(mm)
  }
  out <- eval_("mean", par, mm, control$maxit_l, control$tol_l, control$th, w)
  out$df <- p
  out$p.value <- pchisq(out$statistic, df = out$df, lower.tail = FALSE)
  # out$par <- par
  out$npar <- p
  if (!is.null(weights)) {
    out$weights <- w
  }
  if (model) {
    out$data.matrix <- mm
  }
  out$coefficients <- cf
  class(out) <- "el"
  out
}
