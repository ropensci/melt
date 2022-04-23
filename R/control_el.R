#' Control parameters for computation
#'
#' Specifies details of computation of (constrained) empirical likelihood.
#'
#' @param maxit Maximum number of iterations for constrained minimization of
#'   empirical likelihood. Defaults to \code{200}.
#' @param maxit_l Maximum number of iterations of evaluation of empirical
#'   likelihood. Defaults to \code{50}.
#' @param tol Convergence tolerance for the constrained minimization.
#'   Defaults to \code{1e-06}.
#' @param tol_l Relative convergence tolerance for the evaluation. Defaults
#'   to \code{1e-06}.
#' @param step Step size for projected gradient method. Defaults to \code{NULL}
#'   and set to the reciprocal of sample size.
#' @param th Threshold for negative empirical log-likelihood ratio value.
#'   The iteration stops if the value exceeds the threshold. Defaults to
#'   \code{NULL}.
#' @param nthreads Number of threads for parallel computation via OpenMP (if
#'     available).
#' @details Let \eqn{X_i} be independent and identically distributed
#'   \eqn{p}-dimensional random variables from an unknown distribution \eqn{F}
#'   for \eqn{i = 1, \dots, n}. For a parameter of interest \eqn{\theta(F) \in
#'   {\rm{I\!R}}^p}, consider a \eqn{p}-dimensional smooth estimating function
#'   \eqn{g(X_i, \theta)} with a moment condition
#'   \deqn{\textnormal{E}[g(X_i, \theta)] = 0.}
#'   We assume that there exists an unique \eqn{\theta_0} that solves the above
#'   equation. Given a value of \eqn{\theta}, the (profile) empirical likelihood
#'   ratio is defined by
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
#'   Then the empirical log-likelihood ratio is given by
#'   \deqn{\log\mathcal{R}(\theta) =  \max_{\lambda}\sum_{i = 1}^n
#'   \log(1 + \lambda^\top g(X_i, \theta)).}
#'   Let \eqn{l(\theta)} denote the minus twice empirical log-likelihood ratio
#'   function. Under some regularity conditions, it is known that
#'   \eqn{l(\theta_0)} converges in distribution to \eqn{\chi^2_p}, where
#'   \eqn{\chi^2_p} has a chi-squared distribution with \eqn{p} degrees of
#'   freedom.
#'
#'   Inference for \eqn{\theta} often involves a hypothesis testing
#'   through a constraint on the parameter space. Consider a linear hypothesis
#'   of the form \deqn{L\theta = r,} where the left-hand-side \eqn{L} is a
#'   \eqn{q} by \eqn{p} matrix and the right-hand-side \eqn{r} is a
#'   \eqn{q}-dimensional vector. With additional conditions,
#'   \eqn{l(\theta)} converges in distribution to \eqn{\chi^2_q} under the
#'   constraint of hypothesis, i.e.,
#'   \deqn{\min_{\theta: L\theta = r} l(\theta) \to_d \chi^2_q .}
#'
#'   Most of the functions provided by \strong{melt} perform either of the two
#'   tasks:
#'   \itemize{
#'   \item Evaluating \eqn{l(\theta)} for a specific value of \eqn{\theta},
#'   e.g., \code{\link{el_eval}} and \code{\link{el_mean}}.
#'   \item Minimizing \eqn{l(\theta)} subject to a constraint, e.g.,
#'   \code{\link{el_lm}} and \code{\link{lht}}.
#'   }
#'   As described above, evaluating \eqn{l(\theta)} involves optimization with
#'   respect to \eqn{\lambda}. This problem can be efficiently solved by the
#'   Newton-Raphson method when the interior of convex hull constructed from
#'   the set of \eqn{g(X_i, \theta)} contains the zero vector.
#'
#'   Minimization of \eqn{l(\theta)} with respect to \eqn{\theta}, on the other
#'   hand, is computationally expensive since it implicitly involves the
#'   evaluation step. Further, depending on the form of \eqn{g(X_i, \theta)} and
#'   the constraint, the optimization problem can be nonconvex and have multiple
#'   local minima. For this reason, \strong{melt} only considers linear
#'   hypotheses and performs local minimization of \eqn{l(\theta)} using
#'   projected gradient descent method. With the orthogonal projection matrix
#'   \eqn{P} and a step size \eqn{\gamma}, the algorithm updates \eqn{\theta} as
#'   \deqn{\theta^{(k + 1)} \leftarrow \theta^{(k)} -
#'   \gamma P \nabla l(\theta^{(k)}),}
#'   where \eqn{\nabla l(\theta^{(k)})} denotes the gradient of \eqn{l} at
#'   \eqn{\theta^{(k)}}. The first order optimality condition is
#'   \eqn{P \nabla l(\theta) = 0}, which is used as the stopping criterion.
#' @return A list of class \code{"control_el"} that specifies details of the
#'   optimization with respect to \eqn{\lambda} and \eqn{\theta} with the
#'   following components:
#'   \item{maxit}{Maximum number of iterations for the optimization with
#'   respect to \eqn{\theta}.}
#'   \item{maxit_l}{Maximum number of iterations for the optimization with
#'   respect to \eqn{\lambda}.}
#'   \item{tol}{Convergence tolerance denoted by \eqn{\epsilon}. The iteration
#'   stops when
#'   \deqn{\|P \nabla l(\theta^{(k)})\| < \epsilon.}}
#'   \item{tol_l}{Relative convergence tolerance denoted by \eqn{\delta}. The
#'   iteration stops when
#'   \deqn{\|\lambda^{(k)} - \lambda^{(k - 1)}\| <
#'   \delta\|\lambda^{(k - 1)}\| + \delta^2.}}
#'   \item{step}{Step size \eqn{\gamma} for the projected gradient descent
#'   method.}
#'   \item{th}{Threshold for the negative empirical log-likelihood ratio value.
#'   The iteration stops if the value exceeds the threshold. Defaults to
#'   \code{NULL} and sets the threshold to \code{200 * d}, where \code{d}
#'   corresponds to the degrees of freedom of the limiting chi-squared
#'   distribution of the statistic.}
#'   \item{nthreads}{Number of threads for parallel computation via OpenMP (if
#'     available). Defaults to the half of the available threads. For better
#'     performance, it is recommended to limit the number of threads to the
#'     number of physical cores. Note that it only applies to the following
#'     functions that involve multiple evaluations or minimizations:
#'     \itemize{
#'     \item{\code{\link{confreg}}}
#'     \item{\code{\link{el_lm}}}
#'     \item{\code{\link{el_glm}}}}}
#' @references Adimari, Gianfranco, and Annamaria Guolo. 2010.
#'   “A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.”
#'   Statistical Methods & Applications 19 (4): 463–76.
#'   \doi{10.1007/s10260-010-0137-9}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.”
#'   The Annals of Statistics 19 (4): 1725–47. \doi{10.1214/aos/1176348368}.
#' @references Qin, Jin, and Jerry Lawless. 1994.
#'   “Empirical Likelihood and General Estimating Equations.”
#'   The Annals of Statistics 22 (1): 300–325. \doi{10.1214/aos/1176325370}.
#' @references Qin, Jing, and Jerry Lawless. 1995.
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics 23 (2): 145–59. \doi{10.2307/3315441}.
#' @seealso \link{el_eval}, \link{lht}
#' @examples
#' optcfg <- control_el(maxit = 300L, th = 200, nthreads = 1L)
#' @export
control_el <- function(maxit = 200L, maxit_l = 50L, tol = 1e-06, tol_l = 1e-06,
                       step = NULL, th = NULL, nthreads) {
  # maxit: integer (positive)
  maxit <- tryCatch(as.integer(maxit), warning = function(w) NA,
                    error = function(e) NA)
  if (any(length(maxit) != 1L, is.na(maxit))) {
    stop("'maxit' is not an integer")
  }
  if (maxit < 1) {
    stop("'maxit' is not a positive integer")
  }
  # maxit_l: integer (positive)
  maxit_l <- tryCatch(as.integer(maxit_l), warning = function(w) NA,
                      error = function(e) NA)
  if (any(length(maxit_l) != 1L, is.na(maxit_l))) {
    stop("'maxit_l' is not an integer")
  }
  if (maxit_l < 1) {
    stop("'maxit_l' is not a positive integer")
  }
  # tol: numeric (positive, finite)
  tol <- tryCatch(as.numeric(tol), warning = function(w) NA,
                  error = function(e) NA)
  if (any(length(tol) != 1L, is.na(tol), is.infinite(tol))) {
    stop("'tol' is not a number")
  }
  if (tol < .Machine$double.eps) {
    stop("'tol' is too small")
  }
  # tol_l: numeric (positive, finite)
  tol_l <- tryCatch(as.numeric(tol_l), warning = function(w) NA,
                    error = function(e) NA)
  if (any(length(tol_l) != 1L, is.na(tol_l), is.infinite(tol_l))) {
    stop("'tol' is not a number")
  }
  if (tol_l < .Machine$double.eps) {
    stop("'tol' is too small")
  }
  # step: numeric (positive, finite)
  if (!is.null(step)) {
    step <- tryCatch(as.numeric(step), warning = function(w) NA,
                     error = function(e) NA)
    if (any(length(step) != 1L, is.na(step), is.infinite(step))) {
      stop("'step' is not a number")
    }
    if (step < .Machine$double.eps) {
      stop("'step' is too small")
    }
  }
  # th: numeric (positive, finite)
  if (!is.null(th)) {
    th <- tryCatch(as.numeric(th), warning = function(w) NA,
                   error = function(e) NA)
    if (any(length(th) != 1L, is.na(th), is.infinite(th))) {
      stop("'th' is not a number")
    }
    if (th < .Machine$double.eps) {
      stop("'th' is too small")
    }
  }
  # nthreads: integer (positive)
  max_threads <- max_threads_()
  if (missing(nthreads)) {
    nthreads <- max(1L, max_threads / 2L)
  } else {
    nthreads <- tryCatch(as.integer(nthreads), warning = function(w) NA,
                         error = function(e) NA)
    if (any(length(nthreads) != 1L, is.na(nthreads))) {
      stop("'nthreads' is not an integer")
    }
    if (nthreads < 1) {
      warning("'nthreads' is set to 1")
      nthreads <- 1L
    }
    if (nthreads > max_threads) {
      warning("'nthreads' is set to the maximum number of threads available")
      nthreads <- max_threads
    }
  }
  out <- list(maxit = maxit, maxit_l = maxit_l, tol = tol, tol_l = tol_l,
              step = step, th = th, nthreads = nthreads)
  class(out) <- "control_el"
  out
}
