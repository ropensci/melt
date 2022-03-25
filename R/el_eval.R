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
#' @details Let \eqn{X_i \in {\rm{I\!R}}^p} be i.i.d. random variables for
#'   \eqn{i = 1, \dots, n}. Assume that there exists an unique \eqn{\theta_0 \in
#'   {\rm{I\!R}}^p} that solves \eqn{\textnormal{E}[g(X_i, \theta)] = 0}, where
#'   the estimating function \eqn{g(X_i, \theta)} takes values in
#'   \eqn{{\rm{I\!R}}^p}. Given a value of \eqn{\theta}, the empirical
#'   likelihood ratio is obtained by
#'   \deqn{\mathcal{R}(\theta) =
#'   \max_{p_i}\left\{\prod_{i = 1}^n np_i :
#'   \sum_{i = 1}^n p_i g(X_i, \theta) = 0, p_i \geq 0, \sum_{i = 1}^n p_i = 1
#'   \right\}.}
#'   The Lagrange multiplier \eqn{\lambda \equiv \lambda(\theta)} of the dual
#'   problem leads to
#'   \deqn{p_i = \frac{1}{n}\frac{1}{1 + \lambda^\top g(X_i, \theta)},}
#'   where \eqn{\lambda} solves
#'   \deqn{\frac{1}{n}\sum_{i = 1}^n \frac{g(X_i, \theta)}
#'   {1 + \lambda^\top g(X_i, \theta)} = 0.}
#'   Then the log empirical likelihood ratio is given by
#'   \deqn{\log\mathcal{R}(\theta) = -\sum_{i = 1}^n
#'   \log(1 + \lambda^\top g(X_i, \theta)).}
#'
#'   \code{el_eval} performs the optimization via Newton’s algorithm to compute
#'   \eqn{\lambda} with a \eqn{n} by \eqn{p} numeric matrix argument \code{g},
#'   whose \eqn{i}th is \eqn{g(X_i, \theta)}. If \code{weights} is non
#'   \code{NULL}, the weights are rescaled to add up to \eqn{n}. The
#'   \code{control} argument is a list that can supply any of the following
#'   components:
#' \describe{
#'   \item{maxit}{The maximum number of iterations for the optimization.
#'   Defaults to \code{100}.}
#'   \item{tol}{The relative convergence tolerance, denoted by \eqn{\epsilon}.
#'   The iteration stops when
#'   \deqn{\|\lambda_{k} - \lambda_{k - 1}\| \leq
#'   \epsilon\|\lambda_{k - 1}\| + \epsilon.} Defaults to \code{1e-06}.}
#'   \item{th}{The threshold for the negative log empirical likelihood
#'   ratio value. The iteration stops if the value exceeds the threshold.
#'   Defaults to \code{NULL} and sets the threshold to \eqn{20p}.}
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
#'       \item{convergence}{A logical vector. \code{TRUE} indicates
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
    out <- EL_eval(g, optcfg$maxit, optcfg$tol, optcfg$th)
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
    out <- WEL_eval(g, w, optcfg$maxit, optcfg$tol, optcfg$th)
  }
  out
}
