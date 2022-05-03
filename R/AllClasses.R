#' S4 class \linkS4class{EL}
#'
#' S4 class for empirical likelihood.
#'
#' @details Let \eqn{X_i} be independent and identically distributed
#'   \eqn{p}-dimensional random variable from an unknown distribution \eqn{F}
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
#'   This problem can be efficiently solved by the Newton-Raphson method when
#'   the zero vector is contained in the interior of the convex hull of
#'   \eqn{\{g(X_i, \theta)\}_{i = 1}^n}.
#'
#'   Under some regularity conditions, it is known that
#'   \eqn{-2\log\mathcal{R}(\theta_0)} converges in distribution to
#'   \eqn{\chi^2_p}, where \eqn{\chi^2_p} has a chi-square distribution with
#'   \eqn{p} degrees of freedom.
#' @slot optim A list with the following optimization results:
#'   \itemize{
#'   \item{\code{method } }{Character for method dispatch in internal
#'   functions.}
#'   \item{\code{par } }{Parameter value specified.}
#'   \item{\code{lambda } }{Lagrange multiplier of the dual problem.}
#'   \item{\code{iterations } }{Number of iterations performed.}
#'   \item{\code{convergence } }{Convergence status.}
#'   }
#' @slot logp Log probabilities obtained from empirical likelihood.
#' @slot logl Empirical log-likelihood.
#' @slot loglr Empirical log-likelihood ratio.
#' @slot statistic Minus twice the empirical log-likelihood ratio statistic that
#'   has an asymptotic chi-square distribution.
#' @slot df Degrees of freedom of the statistic.
#' @slot pval \eqn{p}-value of the statistic.
#' @slot npar Number of parameters.
#' @slot weights Rescaled weights used for model fitting.
#' @slot data Data matrix used for model fitting.
#' @slot coefficients Maximum empirical likelihood estimates of the parameters.
#' @references Glenn, N.L., and Yichuan Zhao. 2007.
#'   “Weighted Empirical Likelihood Estimates and Their Robustness Properties.”
#'   Computational Statistics & Data Analysis 51 (10): 5130–41.
#'   \doi{10.1016/j.csda.2006.07.032}.
#' @references Qin, Jin, and Jerry Lawless. 1994.
#'   “Empirical Likelihood and General Estimating Equations.”
#'   The Annals of Statistics 22 (1): 300–325. \doi{10.1214/aos/1176325370}.
#' @examples
#' showClass("EL")
setClass("EL",
  slots = c(
    optim = "list", logp = "numeric", logl = "numeric", loglr = "numeric",
    statistic = "numeric", df = "integer", pval = "numeric", npar = "integer",
    weights = "numeric", data = "matrix", coefficients = "numeric"
  ),
  prototype = list(
    optim = list(), logp = numeric(), logl = numeric(), loglr = numeric(),
    statistic = numeric(), df = 0L, pval = numeric(), npar = 0L,
    weights = numeric(), data = matrix(NA_real_, nrow = 0L, ncol = 0L),
    coefficients = numeric()
  )
  # slots = c(
  #   optim = "list", logp = "numeric", logl = "numeric", loglr = "numeric",
  #   statistic = "numeric", df = "integer", pval = "numeric", npar = "integer",
  #   weights = "numeric", data = "ANY", coefficients = "numeric"
  # ),
  # prototype = list(
  #   optim = list(), logp = numeric(), logl = numeric(), loglr = numeric(),
  #   statistic = numeric(), df = 0L, pval = numeric(), npar = 0L,
  #   weights = numeric(),
  #   coefficients = numeric()
  # )
)

#' S4 class \linkS4class{CEL}
#'
#' S4 class for constrained empirical likelihood. It inherits from
#' \linkS4class{EL} class. Note that \code{optim} slot has constrained
#' optimization results with respect to parameters, not the Lagrange multiplier.
#'
#' @details Let \eqn{l(\theta)} denote the minus twice the empirical
#'   log-likelihood ratio function. We consider a linear hypothesis of the form
#'   \deqn{L\theta = r,} where the left-hand-side \eqn{L} is a \eqn{q} by
#'   \eqn{p} matrix and the right-hand-side \eqn{r} is a \eqn{q}-dimensional
#'   vector. Under some regularity conditions, \eqn{l(\theta)} converges in
#'   distribution to \eqn{\chi^2_q} under the constraint of hypothesis, i.e.,
#'   \deqn{\min_{\theta: L\theta = r} l(\theta) \to_d \chi^2_q .}
#'
#'   Minimization of \eqn{l(\theta)} with respect to \eqn{\theta} is
#'   computationally expensive since it implicitly involves the
#'   evaluation step as described in \linkS4class{EL}. Further, depending on the
#'   form of \eqn{g(X_i, \theta)} and the constraint, the optimization problem
#'   can be nonconvex and have multiple local minima. For this reason,
#'   \strong{melt} only considers linear hypotheses and performs local
#'   minimization of \eqn{l(\theta)} using projected gradient descent method.
#'   With the orthogonal projection matrix \eqn{P} and a step size \eqn{\gamma},
#'   the algorithm updates \eqn{\theta} as
#'   \deqn{\theta^{(k + 1)} \leftarrow \theta^{(k)} -
#'   \gamma P \nabla l(\theta^{(k)}),}
#'   where \eqn{\nabla l(\theta^{(k)})} denotes the gradient of \eqn{l} at
#'   \eqn{\theta^{(k)}}. The first order optimality condition is
#'   \eqn{P \nabla l(\theta) = 0}, which is used as the stopping criterion.
#' @slot optim A list with the following optimization results:
#'   \itemize{
#'   \item{\code{method } }{Character for method dispatch in internal
#'   functions.}
#'   \item{\code{par } }{Parameter value that minimizes the constrained
#'   empirical likelihood.}
#'   \item{\code{lambda } }{Lagrange multiplier of the dual problem.}
#'   \item{\code{iterations } }{Number of iterations performed.}
#'   \item{\code{convergence } }{Convergence status.}
#'   }
#' @slot logp Log probabilities obtained from constrained empirical likelihood.
#' @slot logl Constrained empirical log-likelihood.
#' @slot loglr Constrained empirical log-likelihood ratio.
#' @slot statistic Minus twice the constrained empirical log-likelihood ratio
#' statistic that has an asymptotic chi-square distribution.
#' @slot df Degrees of freedom of the statistic.
#' @slot pval \eqn{p}-value of the statistic.
#' @slot npar Number of parameters.
#' @slot weights Rescaled weights used for model fitting.
#' @slot data Data matrix used for model fitting.
#' @slot coefficients Maximum empirical likelihood estimates of the parameters.
#' @references Adimari, Gianfranco, and Annamaria Guolo. 2010.
#'   “A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.”
#'   Statistical Methods & Applications 19 (4): 463–76.
#'   \doi{10.1007/s10260-010-0137-9}.
#' @references Qin, Jing, and Jerry Lawless. 1995.
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics 23 (2): 145–59. \doi{10.2307/3315441}.
#' @examples
#' showClass("CEL")
setClass("CEL", contains = "EL")

#' S4 class \linkS4class{LM}
#'
#' S4 class for linear models with empirical likelihood. It inherits from
#' \linkS4class{CEL} class.
#'
#' @details If there is no intercept in a model, \code{optim}
#' slot need to be understood in terms of \linkS4class{EL} class since
#' constrained optimization is not involved in the overall test.
#' @slot parTests A list with the test results for each parameter:
#'   \itemize{
#'   \item{\code{statistic } }{A numeric vector of chi-squared statistics.}
#'   \item{\code{convergence } }{Convergence status of tests for each parameter.
#'   }
#'   }
#' @slot misc A list with miscellaneous outputs from a model fitting function.
#' They are used in other generics and methods.
#' @examples
#' showClass("LM")
setClass("LM", contains = "CEL", slots = c(parTests = "list", misc = "list"))

#' S4 class \linkS4class{GLM}
#'
#' S4 class for generalized linear models with empirical likelihood. It inherits
#' from \linkS4class{LM} class.
#'
#' @examples
#' showClass("GLM")
setClass("GLM", contains = "LM")

#' S4 class \linkS4class{ConfregEL}
#'
#' S4 class for confidence region.
#'
#' @slot points A numeric matrix with two columns for boundary points of a
#'   confidence region.
#' @slot estimates A numeric vector of length two for parameter estimates.
#' @slot level A confidence level required.
#' @slot cv A critical value for calibration of empirical likelihood ratio
#'   statistic.
#' @slot pnames A character vector of length two for the name of parameters.
#' @examples
#' showClass("ConfregEL")
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

#' S4 class \linkS4class{ELD}
#'
#' S4 class for empirical likelihood displacement.
#'
#' @slot eld A numeric vector of empirical likelihood displacement values.
#' @examples
#' showClass("ELD")
setClass("ELD", slots = c(eld = "numeric"))

#' S4 class \linkS4class{SummaryLM}
#'
#' S4 class for a summary of \linkS4class{LM} objects.
#'
#' @slot statistic Minus twice the constrained empirical log-likelihood ratio
#'   for the overall test of the model.
#' @slot df Degrees of freedom of the statistic.
#' @slot convergence Convergence status of the minimization.
#' @slot parMatrix A numeric matrix of the test results of the parameters.
#' @slot weighted A logical for whether the given model is weighted or not.
#' @slot na.action Information returned by \code{\link[stats]{model.frame}} on
#'   the special handling of NAs.
#' @slot call Matched call.
#' @slot terms \code{\link[stats]{terms}} object used.
#' @slot aliased Named logical vector showing if the original coefficients are
#'   aliased.
#' @examples
#' showClass("SummaryLM")
setClass("SummaryLM", slots = c(
  statistic = "numeric", df = "integer", convergence = "logical",
  parMatrix = "matrix", weighted = "logical", na.action = "ANY", call = "ANY",
  terms = "ANY", aliased = "logical"
))

#' S4 class \linkS4class{ControlEL}
#'
#' S4 class for details of computation of empirical likelihood.
#'
#' @slot maxit Maximum number of iterations for the optimization with
#'   respect to \eqn{\theta}.
#' @slot maxit_l Maximum number of iterations for the optimization with
#'   respect to \eqn{\lambda}.
#' @slot tol Convergence tolerance denoted by \eqn{\epsilon}. The iteration
#'   stops when
#'   \deqn{\|P \nabla l(\theta^{(k)})\| < \epsilon.}
#' @slot tol_l Relative convergence tolerance denoted by \eqn{\delta}. The
#'   iteration stops when
#'   \deqn{\|\lambda^{(k)} - \lambda^{(k - 1)}\| <
#'   \delta\|\lambda^{(k - 1)}\| + \delta^2.}
#' @slot step Step size \eqn{\gamma} for the projected gradient descent
#'   method.
#' @slot th Threshold for the negative empirical log-likelihood ratio value.
#'   The iteration stops if the value exceeds the threshold. Defaults to
#'   \code{NULL} and sets the threshold to \code{200 * d}, where \code{d}
#'   corresponds to the degrees of freedom of the limiting chi-squared
#'   distribution of the statistic.
#' @slot nthreads Number of threads for parallel computation via OpenMP (if
#'   available). Defaults to the half of the available threads. For better
#'   performance, it is recommended to limit the number of threads to the
#'   number of physical cores. Note that it only applies to the following
#'   functions that involve multiple evaluations or minimizations:
#'   \itemize{
#'   \item{\code{\link{confreg}}}
#'   \item{\code{\link{el_lm}}}
#'   \item{\code{\link{el_glm}}}
#'   \item{\code{\link{eld}}}}
#' @examples
#' showClass("ControlEL")
setClass("ControlEL",
  slots = c(
    maxit = "integer", maxit_l = "integer", tol = "numeric", tol_l = "numeric",
    step = "ANY", th = "ANY", nthreads = "integer"
  ),
  prototype = list(
    maxit = 200L, maxit_l = 50L, tol = 1e-06, tol_l = 1e-06,
    step = NULL, th = NULL, nthreads = NULL
  )
)

#' S4 class \linkS4class{logLikEL}
#'
#' S4 class for empirical log-likelihood.
#'
#' @slot logLik Empirical log-likelihood.
#' @slot df Degrees of freedom or the number of (estimated) parameters in the
#'   model.
#' @examples
#' showClass("logLikEL")
setClass("logLikEL", slots = c(logLik = "numeric", df = "integer"))
