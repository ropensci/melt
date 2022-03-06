#' Empirical likelihood with general estimating functions
#'
#' Computes empirical likelihood with general estimating functions.
#'
#' @param g A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation of an estimating
#'   function.
#' @param weights An optional numeric vector of weights.
#'   Defaults to \code{NULL}, corresponding to identical weights.
#'   If non \code{NULL}, weighted empirical likelihood is computed.
#' @param control A list of control parameters. See ‘Details’.
#' @details For a \eqn{p}-dimensional parameter \eqn{\theta_0} and an estimating
#'   function \eqn{g}, consider the following moment condition:
#'   \deqn{E{g(X_i, \theta_0)} = 0, i = 1, \dots, n.}
#'   Given a value of \eqn{\theta}, \code{el_eval} computes the empirical
#'   likelihood with a \eqn{n} by \eqn{p} numeric matrix argument \code{g},
#'   whose rows consist of all \eqn{g(X_i, \theta)}. If \code{weights} is non
#'   \code{NULL}, the weights are rescaled to add up to \eqn{n}. The
#'   \code{control} argument is a list that can supply any of the following
#'   components:
#' \describe{
#'   \item{maxit}{The maximum number of iterations for optimization. Defaults to
#'   \code{100}.}
#'   \item{abstol}{The absolute convergence tolerance for log likelihood ratio
#'   value. Defaults to \code{1e-06}.}
#'   \item{threshold}{The threshold for log likelihood ratio value. The
#'   computation stops if the value exceeds the threshold. Defaults to
#'   \code{NULL} and sets the threshold to \eqn{20p}.}
#' }
#' @return A list with the following components:
#' \describe{
#'   \item{optim}{A list with the following optimization results:
#'     \describe{
#'       \item{lambda}{The Lagrange multiplier of dual problem.}
#'       \item{weights}{If non \code{NULL} \code{weights} is supplied, the
#'       rescaled weights are returned.}
#'       \item{logLR}{The (weighted) log empirical likelihood ratio value.}
#'       \item{iterations}{The number of iterations performed.}
#'       \item{convergence}{A logical vector. \code{TRUE} indicates the
#'       convergence of the algorithm.}
#'     }
#'   }
#'   \item{statistic}{The chi-square statistic.}
#'   \item{df}{The degrees of freedom of the statistic.}
#'   \item{p.value}{The \eqn{p}-value of the statistic.}
#' }
#' @references Glenn, N.L., and Yichuan Zhao. 2007.
#'   “Weighted Empirical Likelihood Estimates and Their Robustness Properties.”
#'   Computational Statistics & Data Analysis 51 (10): 5130–41.
#'   \doi{10.1016/j.csda.2006.07.032}.
#' @references Qin, Jin, and Jerry Lawless. 1994.
#'   “Empirical Likelihood and General Estimating Equations.”
#'   The Annals of Statistics 22 (1).
#'   \doi{10.1214/aos/1176325370}.
#' @export
el_eval <- function(g, weights = NULL, control = list())
{
  # check g
  if (!is.matrix(g))
    g <- as.matrix(g)
  if (!is.numeric(g) || any(!is.finite(g)))
    stop("'g' must be a finite numeric matrix")
  if (NROW(g) < 2L)
    stop("not enough 'g' observations")

  # check control
  optcfg <- check_control(control)
  if (is.null(weights)) {
    out <- EL_eval(g, optcfg$maxit, optcfg$abstol, optcfg$threshold)
  } else {
    if (!is.numeric(weights))
      stop("'weights' must be a numeric vector")
    w <- as.numeric(weights)
    if (any(!is.finite(w)))
      stop("'weights' must be a finite numeric vector")
    if (any(w < 0))
      stop("negative 'weights' are not allowed")
    if (length(w) != NROW(g))
      stop("'g' and 'weights' have incompatible dimensions")
    w <- (NROW(g) / sum(w)) * w
    out <- WEL_eval(g, w, optcfg$maxit, optcfg$abstol, optcfg$threshold)
  }
  out
}
