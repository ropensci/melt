#' \linkS4class{EL} class
#'
#' S4 class for empirical likelihood.
#'
#' @details Let \eqn{X_i} be independent and identically distributed
#'   \eqn{p}-dimensional random variable from an unknown distribution \eqn{F}
#'   for \eqn{i = 1, \dots, n}. We assume that \eqn{F} has a positive definite
#'   covariance matrix. For a parameter of interest
#'   \eqn{\theta(F) \in {\rm{I\!R}}^p}, consider a \eqn{p}-dimensional smooth
#'   estimating function \eqn{g(X_i, \theta)} with a moment condition
#'   \deqn{\textrm{E}[g(X_i, \theta)] = 0.}
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
#'   \deqn{\log\mathcal{R}(\theta) = -\sum_{i = 1}^n
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
#'   \item{`par ` }{A numeric vector of the specified parameters.}
#'   \item{`lambda ` }{A numeric vector of the Lagrange multipliers.}
#'   \item{`iterations ` }{A single integer for the number of iterations
#'   performed.}
#'   \item{`convergence ` }{A single logical for the convergence status.}
#'   }
#' @slot logp A numeric vector of the log probabilities obtained from empirical
#'   likelihood.
#' @slot logl A single numeric for the empirical log-likelihood.
#' @slot loglr A single numeric for the empirical log-likelihood ratio.
#' @slot statistic A single numeric for the minus twice the empirical
#'   log-likelihood ratio statistic that has an asymptotic chi-square
#'   distribution.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot npar A single integer for the number of parameters.
#' @slot weights A numeric vector of re-scaled weights used for model fitting.
#' @slot data A numeric matrix for the data used for model fitting.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @aliases EL
#' @references Owen A (2001). Empirical Likelihood. Chapman & Hall/CRC.
#'   \doi{10.1201/9781420036152}.
#' @references Qin J, Lawless J (1994).
#'   “Empirical Likelihood and General Estimating Equations.”
#'   The Annals of Statistics, 22(1), 300–325. \doi{10.1214/aos/1176325370}.
#' @examples
#' showClass("EL")
#' @srrstats {G1.0} References are included throughout the package
#'   documentation.
#' @srrstats {G1.4} `roxygen2` is used to document all functions.
#' @srrstats {G1.3} Statistical terminology is unambiguously defined throughout
#'   the package documentation.
#' @srrstats {RE1.4} Some of the core distributional assumptions are documented.
#'   We note that these are not exhaustive and added the phrase "Under some
#'   regularity conditions" instead. For more details, interested users are
#'   expected to resort to published references. Additionally, violation of the
#'   convex hull constraint is tested in `test/testthat/test-el_mean.R`.
setClass("EL",
  slots = c(
    optim = "list", logp = "numeric", logl = "numeric", loglr = "numeric",
    statistic = "numeric", df = "integer", pval = "numeric", nobs = "integer",
    npar = "integer", weights = "numeric", data = "ANY",
    coefficients = "numeric", method = "character"
  ),
  prototype = list(
    optim = list(), logp = numeric(), logl = numeric(), loglr = numeric(),
    statistic = numeric(), df = 0L, pval = numeric(), npar = 0L,
    weights = numeric(), data = NULL,
    coefficients = numeric(), method = NA_character_
  )
)


#' \linkS4class{CEL} class
#'
#' S4 class for constrained empirical likelihood. It inherits from
#' \linkS4class{EL} class. Note that `optim` slot has constrained optimization
#' results with respect to parameters, not the Lagrange multiplier.
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
#'   can be nonconvex and have multiple local minima. For this reason, the
#'   package \pkg{melt} only considers linear hypotheses and performs local
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
#'   \item{`par ` }{A numeric vector of the parameter value that minimizes the
#'   empirical likelihood subject to the constraints.}
#'   \item{`lambda ` }{A numeric vector of the Lagrange multipliers.}
#'   \item{`iterations ` }{A single integer for the number of iterations
#'   performed.}
#'   \item{`convergence ` }{A single logical for the convergence status.}
#'   }
#' @aliases CEL
#' @references Adimari G, Guolo A (2010).
#'   “A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.”
#'   Statistical Methods & Applications, 19(4), 463–476.
#'   \doi{10.1007/s10260-010-0137-9}.
#' @references Qin J, Lawless J (1995).
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics, 23(2), 145–159. \doi{10.2307/3315441}.
#' @examples
#' showClass("CEL")
#' @srrstats {G1.1} The package attempts the first implementation of the nested
#'   bilevel optimization approach within R to compute constrained empirical
#'   likelihood. The inner layer Newton-Raphson method for empirical likelihood
#'   is written in C++, enabling faster computation than other routines written
#'   in R.
setClass("CEL", contains = "EL")


#' \linkS4class{LM} class
#'
#' S4 class for linear models with empirical likelihood. It inherits from
#' \linkS4class{CEL} class.
#'
#' @details If there is no intercept in a model, `optim` slot need to be
#'   understood in terms of \linkS4class{EL} class since constrained
#'   optimization is not involved in the overall test.
#' @slot parTests A list with the test results for each parameter:
#'   \itemize{
#'   \item{`statistic ` }{A numeric vector of the empirical likelihood
#'   ratio statistics.}
#'   \item{`convergence ` }{A logical vector of the convergence status of
#'   tests for each parameter.}
#'   }
#' @slot misc A list with miscellaneous outputs from a model fitting function.
#'   They are used in other generics and methods.
#' @aliases LM
#' @examples
#' showClass("LM")
setClass("LM",
  contains = "CEL",
  slots = c(parTests = "list", call = "call", terms = "ANY", misc = "list")
)

#' \linkS4class{GLM} class
#'
#' S4 class for generalized linear models with empirical likelihood. It inherits
#' from \linkS4class{LM} class.
#'
#' @aliases GLM
#' @examples
#' showClass("GLM")
setClass("GLM", contains = "LM")


#' \linkS4class{SD} class
#'
#' S4 class for standard deviation. It inherits from \linkS4class{EL} class.
#'
#' @aliases SD
#' @examples
#' showClass("SD")
setClass("SD", contains = "EL")


#' \linkS4class{ConfregEL} class
#'
#' S4 class for confidence region.
#'
#' @slot points A numeric matrix with two columns for boundary points of a
#'   confidence region.
#' @slot estimates A numeric vector of length two for parameter estimates.
#' @slot level A single numeric for the confidence level required.
#' @slot cv A single numeric for the critical value for calibration of empirical
#'   likelihood ratio statistic.
#' @slot pnames A character vector of length two for the name of parameters.
#' @aliases ConfregEL
#' @examples
#' showClass("ConfregEL")
setClass("ConfregEL",
  slots = c(
    points = "matrix", estimates = "numeric", level = "numeric", cv = "numeric",
    pnames = "character"
  ),
  prototype = list(
    points = NULL, estimates = NA_real_, cv = NA_real_, pnames = NA_character_
  )
)


#' \linkS4class{ControlEL} class
#'
#' S4 class for computational details of empirical likelihood.
#'
#' @slot maxit A single integer for the maximum number of iterations for the
#'   optimization with respect to \eqn{\theta}.
#' @slot maxit_l A single integer for the maximum number of iterations for the
#'   optimization with respect to \eqn{\lambda}.
#' @slot tol A single numeric for the convergence tolerance denoted by
#'   \eqn{\epsilon}. The iteration stops when
#'   \deqn{\|P \nabla l(\theta^{(k)})\| < \epsilon.}
#' @slot tol_l A single numeric for the relative convergence tolerance denoted
#'   by \eqn{\delta}. The iteration stops when
#'   \deqn{\|\lambda^{(k)} - \lambda^{(k - 1)}\| <
#'   \delta\|\lambda^{(k - 1)}\| + \delta^2.}
#' @slot step A single numeric for the step size \eqn{\gamma} for the projected
#'   gradient descent method.
#' @slot th A single numeric for the threshold for the negative empirical
#'   log-likelihood ratio.
#' @slot verbose A single logical for whether to print a message on the
#'   convergence status.
#' @slot keep_data A single logical for whether to
#' @slot nthreads A single integer for the number of threads for parallel
#'   computation via OpenMP (if available).
#' @slot seed A single integer for the seed for random number generation.
#' @slot b A single integer for the number of bootstrap replicates.
#' @slot m A single integer for the number of Monte Carlo samples.
#' @aliases ControlEL
#' @seealso [el_control()]
#' @examples
#' showClass("ControlEL")
setClass("ControlEL",
  slots = c(
    maxit = "integer", maxit_l = "integer", tol = "numeric", tol_l = "numeric",
    step = "ANY", th = "ANY", verbose = "logical", keep_data = "logical",
    nthreads = "integer", seed = "integer", b = "integer", m = "integer"
  )
)


#' \linkS4class{ELD} class
#'
#' S4 class for empirical likelihood displacement.
#'
#' @slot eld A numeric vector of empirical likelihood displacement values.
#' @aliases ELD
#' @examples
#' showClass("ELD")
setClass("ELD", slots = c(eld = "numeric"))


#' \linkS4class{ELT} class
#'
#' S4 class for empirical likelihood test.
#'
#' @slot optim A list with the optimization results.
#' @slot alpha A single numeric for the significance level.
#' @slot logl A single numeric for the (constrained) empirical log-likelihood.
#' @slot loglr A single numeric for the (constrained) empirical log-likelihood
#'   ratio.
#' @slot statistic A single numeric for the minus twice the (constrained)
#'   empirical log-likelihood ratio.
#' @slot cv A single numeric for the critical value.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot calibrate A single character for the calibration method used.
#' @aliases ELT
#' @examples
#' showClass("ELT")
setClass("ELT",
  slots = c(
    optim = "list", alpha = "numeric", logl = "numeric", loglr = "numeric",
    statistic = "numeric", cv = "numeric", pval = "numeric",
    calibrate = "character"
  )
)


#' \linkS4class{ELMT} class
#'
#' S4 class for empirical likelihood multiple tests.
#'
#' @slot alpha A single numeric for the overall significance level.
#' @slot statistic A numeric vector for the minus twice the (constrained)
#'   empirical log-likelihood ratios.
#' @slot cv A single numeric for the multiplicity adjusted critical value.
#' @slot pval A numeric vector for the multiplicity adjusted \eqn{p}-values.
#' @slot calibrate A single character for the calibration method used.
#' @aliases ELMT
#' @examples
#' showClass("ELMT")
setClass("ELMT",
  slots = c(
    alpha = "numeric", statistic = "numeric", cv = "numeric", pval = "numeric",
    calibrate = "character"
  )
)


#' \linkS4class{logLikEL} class
#'
#' S4 class for empirical log-likelihood.
#'
#' @slot logLik A single numeric for the empirical log-likelihood.
#' @slot df A single integer for the degrees of freedom or the number of
#'   (estimated) parameters in the model.
#' @aliases logLikEL
#' @examples
#' showClass("logLikEL")
setClass("logLikEL", slots = c(logLik = "numeric", df = "integer"))


#' \linkS4class{SummaryLM} class
#'
#' S4 class for a summary of \linkS4class{LM} objects.
#'
#' @slot statistic A single numeric for the minus twice the empirical
#'   log-likelihood ratio for the overall test of the model.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot convergence A single logical for the convergence status of the
#'   constrained minimization.
#' @slot parMatrix A numeric matrix of the test results of the parameters.
#' @slot weighted A single logical for whether the given model is weighted or
#'   not.
#' @slot na.action Information returned by [`model.frame`] on the special
#'   handling of `NA`s.
#' @slot call Matched call.
#' @slot terms [`terms`] object used.
#' @slot aliased A named logical vector showing if the original coefficients are
#'   aliased.
#' @aliases SummaryLM
#' @examples
#' showClass("SummaryLM")
setClass("SummaryLM", slots = c(
  statistic = "numeric", df = "integer", convergence = "logical",
  parMatrix = "matrix", weighted = "logical", na.action = "ANY", call = "ANY",
  terms = "ANY", aliased = "logical"
))
