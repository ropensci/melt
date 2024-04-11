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
#' @slot keep_data A single logical for whether to keep the data used for
#'   fitting model objects.
#' @slot nthreads A single integer for the number of threads for parallel
#'   computation via OpenMP (if available).
#' @slot seed A single integer for the seed for random number generation.
#' @slot an A single numeric representing the scaling factor for adjusted
#'   empirical likelihood calibration.
#' @slot b A single integer for the number of bootstrap replicates.
#' @slot m A single integer for the number of Monte Carlo samples.
#' @aliases ControlEL
#' @examples
#' showClass("ControlEL")
setClass("ControlEL",
  slots = c(
    maxit = "integer", maxit_l = "integer", tol = "numeric", tol_l = "numeric",
    step = "ANY", th = "ANY", verbose = "logical", keep_data = "logical",
    nthreads = "integer", seed = "ANY", an = "ANY", b = "integer", m = "integer"
  )
)

#' \linkS4class{EL} class
#'
#' S4 class for empirical likelihood.
#'
#' @details Let \eqn{X_i} be independent and identically distributed
#'   \eqn{p}-dimensional random variable from an unknown distribution \eqn{P}
#'   for \eqn{i = 1, \dots, n}. We assume that \eqn{P} has a positive definite
#'   covariance matrix. For a parameter of interest
#'   \eqn{\theta(F) \in {\rm{I\!R}}^p}, consider a \eqn{p}-dimensional smooth
#'   estimating function \eqn{g(X_i, \theta)} with a moment condition
#'   \deqn{\textrm{E}[g(X_i, \theta)] = 0.}
#'   We assume that there exists an unique \eqn{\theta_0} that solves the above
#'   equation. Given a value of \eqn{\theta}, the (profile) empirical likelihood
#'   ratio is defined by
#'   \deqn{R(\theta) =
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
#'   \deqn{\log R(\theta) = -\sum_{i = 1}^n
#'   \log(1 + \lambda^\top g(X_i, \theta)).}
#'   This problem can be efficiently solved by the Newton-Raphson method when
#'   the zero vector is contained in the interior of the convex hull of
#'   \eqn{\{g(X_i, \theta)\}_{i = 1}^n}.
#'
#'   It is known that \eqn{-2\log R(\theta_0)} converges in
#'   distribution to \eqn{\chi^2_p}, where \eqn{\chi^2_p} has a chi-square
#'   distribution with \eqn{p} degrees of freedom. See the references below for
#'   more details.
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the specified parameters.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logp A numeric vector of the log probabilities of the empirical
#'   likelihood.
#' @slot logl A single numeric of the empirical log-likelihood.
#' @slot loglr A single numeric of the empirical log-likelihood ratio.
#' @slot statistic A single numeric of minus twice the empirical log-likelihood
#'   ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weights A numeric vector of the re-scaled weights used for the model
#'   fitting.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot data A numeric matrix of the data for the model fitting.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases EL
#' @references
#'   Owen A (2001).
#'   \emph{Empirical Likelihood}. Chapman & Hall/CRC.
#'   \doi{10.1201/9781420036152}.
#' @references
#'   Qin J, Lawless J (1994).
#'   ``Empirical Likelihood and General Estimating Equations.''
#'   \emph{The Annals of Statistics}, **22**(1), 300--325.
#'   \doi{10.1214/aos/1176325370}.
#' @examples
#' showClass("EL")
setClass("EL",
  slots = c(
    optim = "list", logp = "numeric", logl = "numeric", loglr = "numeric",
    statistic = "numeric", df = "integer", pval = "numeric", nobs = "integer",
    npar = "integer", weights = "numeric", coefficients = "numeric",
    method = "character", data = "ANY", control = "ControlEL"
  )
)


#' \linkS4class{CEL} class
#'
#' S4 class for constrained empirical likelihood. It inherits from
#'   \linkS4class{EL} class. Note that the `optim` slot has constrained
#'   optimization results with respect to the parameters, not the Lagrange
#'   multiplier.
#'
#' @details Let \eqn{l(\theta)} denote minus twice the empirical log-likelihood
#'   ratio function. We consider a linear hypothesis of the form
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
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the constrained optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logp A numeric vector of the log probabilities of the constrained
#'   empirical likelihood.
#' @slot logl A single numeric of the constrained empirical log-likelihood.
#' @slot loglr A single numeric of the constrained empirical log-likelihood
#'   ratio.
#' @slot statistic A single numeric of minus twice the constrained empirical
#'   log-likelihood ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weights A numeric vector of the re-scaled weights used for the model
#'   fitting.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot data A numeric matrix of the data for the model fitting.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases CEL
#' @references
#'   Adimari G, Guolo A (2010).
#'   ``A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.''
#'   \emph{Statistical Methods & Applications}, **19**(4), 463--476.
#'   \doi{10.1007/s10260-010-0137-9}.
#' @references
#'   Qin J, Lawless J (1995).
#'   ``Estimating Equations, Empirical Likelihood and Constraints on
#'   Parameters.'' \emph{Canadian Journal of Statistics}, **23**(2), 145--159.
#'   \doi{10.2307/3315441}.
#' @examples
#' showClass("CEL")
setClass("CEL", contains = "EL")


setOldClass("terms")
#' \linkS4class{LM} class
#'
#' S4 class for linear models with empirical likelihood. It inherits from
#'   \linkS4class{CEL} class.
#'
#' @details The overall test involves a constrained optimization problem. All
#'   the parameters except for the intercept are constrained to zero. The
#'   `optim` slot contains the results. When there is no intercept, all
#'   parameters are set to zero, and the results need to be understood in terms
#'   of \linkS4class{EL} class since no constrained optimization is involved.
#'   Once the solution is found, the log probabilities (`logp`) and the
#'   (constrained) empirical likelihood values (`logl`, `loglr`, `statistic`)
#'   readily follow, along with the degrees of freedom (`df`) and the
#'   \eqn{p}-value (`pval`). The significance tests for each parameter also
#'   involve constrained optimization problems where only one parameter is
#'   constrained to zero. The `sigTests` slot contains the results.
#' @slot sigTests A list of the following results of significance tests:
#'   * `statistic` A numeric vector of minus twice the (constrained) empirical
#'   log-likelihood ratios with asymptotic chi-square distributions.
#'   * `iterations` An integer vector for the number of iterations performed for
#'   each parameter.
#'   * `convergence` A logical vector for the convergence status of each
#'   parameter.
#' @slot call A matched call.
#' @slot terms A [`terms`] object used.
#' @slot misc A list of various outputs obtained from the model fitting process.
#'   They are used in other generics and methods.
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the (constrained) optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#' @slot logp A numeric vector of the log probabilities of the (constrained)
#'   empirical likelihood.
#' @slot logl A single numeric of the (constrained) empirical log-likelihood.
#' @slot loglr A single numeric of the (constrained) empirical log-likelihood
#'   ratio.
#' @slot statistic A single numeric of minus twice the (constrained) empirical
#'   log-likelihood ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weights A numeric vector of the re-scaled weights used for the model
#'   fitting.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot data A numeric matrix of the data for the model fitting.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases LM
#' @examples
#' showClass("LM")
setClass("LM",
  slots = c(sigTests = "ANY", call = "call", terms = "terms", misc = "list"),
  contains = "CEL"
)


setOldClass("family")
#' \linkS4class{GLM} class
#'
#' S4 class for generalized linear models. It inherits from \linkS4class{LM}
#'   class.
#'
#' @details The overall test involves a constrained optimization problem. All
#'   the parameters except for the intercept are constrained to zero. The
#'   `optim` slot contains the results. When there is no intercept, all
#'   parameters are set to zero, and the results need to be understood in terms
#'   of \linkS4class{EL} class since no constrained optimization is involved.
#'   Once the solution is found, the log probabilities (`logp`) and the
#'   (constrained) empirical likelihood values (`logl`, `loglr`, `statistic`)
#'   readily follow, along with the degrees of freedom (`df`) and the
#'   \eqn{p}-value (`pval`). The significance tests for each parameter also
#'   involve constrained optimization problems where only one parameter is
#'   constrained to zero. The `sigTests` slot contains the results.
#' @slot family A [`family`] object used.
#' @slot dispersion A single numeric for the estimated dispersion parameter.
#' @slot sigTests A list of the following results of significance tests:
#'   * `statistic` A numeric vector of minus twice the (constrained) empirical
#'   log-likelihood ratios with asymptotic chi-square distributions.
#'   * `iterations` An integer vector for the number of iterations performed for
#'   each parameter.
#'   * `convergence` A logical vector for the convergence status of each
#'   parameter.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot call A matched call.
#' @slot terms A [`terms`] object used.
#' @slot misc A list of various outputs obtained from the model fitting process.
#'   They are used in other generics and methods.
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the (constrained) optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#' @slot logp A numeric vector of the log probabilities of the (constrained)
#'   empirical likelihood.
#' @slot logl A single numeric of the (constrained) empirical log-likelihood.
#' @slot loglr A single numeric of the (constrained) empirical log-likelihood
#'   ratio.
#' @slot statistic A single numeric of minus twice the (constrained) empirical
#'   log-likelihood ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weights A numeric vector of the re-scaled weights used for the model
#'   fitting.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot data A numeric matrix of the data for the model fitting.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases GLM
#' @examples
#' showClass("GLM")
setClass("GLM",
  slots = c(family = "family", dispersion = "numeric"),
  contains = "LM"
)


#' \linkS4class{ConfregEL} class
#'
#' S4 class for confidence region. It inherits from `"matrix"`.
#'
#' @slot estimates A numeric vector of length two for the parameter estimates.
#' @slot level A single numeric for the confidence level required.
#' @slot cv A single numeric for the critical value for calibration of empirical
#'   likelihood ratio statistic.
#' @slot pnames A character vector of length two for the name of parameters.
#' @aliases ConfregEL
#' @examples
#' showClass("ConfregEL")
setClass("ConfregEL",
  slots = c(
    estimates = "numeric", level = "numeric", cv = "numeric",
    pnames = "character"
  ),
  contains = "matrix"
)


#' \linkS4class{ELD} class
#'
#' S4 class for empirical likelihood displacement. It inherits from `"numeric"`.
#'
#' @aliases ELD
#' @examples
#' showClass("ELD")
setClass("ELD", contains = "numeric")


#' \linkS4class{ELMT} class
#'
#' S4 class for empirical likelihood multiple tests.
#'
#' @slot estimates A list of numeric vectors of the estimates of the linear
#'   hypotheses.
#' @slot statistic A numeric vector of minus twice the (constrained) empirical
#'   log-likelihood ratios with asymptotic chi-square distributions.
#' @slot df An integer vector of the marginal degrees of freedom of the
#'   statistic.
#' @slot pval A numeric vector for the multiplicity adjusted \eqn{p}-values.
#' @slot cv A single numeric for the multiplicity adjusted critical value.
#' @slot rhs A numeric vector for the right-hand sides of the hypotheses.
#' @slot lhs A numeric matrix for the left-hand side of the hypotheses.
#' @slot alpha A single numeric for the overall significance level.
#' @slot calibrate A single character for the calibration method used.
#' @slot weights A numeric vector of the re-scaled weights used for the model
#'   fitting.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot data A numeric matrix of the data for the model fitting.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases ELMT
#' @examples
#' showClass("ELMT")
setClass("ELMT",
  slots = c(
    estimates = "list", statistic = "numeric", df = "integer", pval = "numeric",
    cv = "numeric", rhs = "numeric", lhs = "matrix", alpha = "numeric",
    calibrate = "character", weights = "numeric", coefficients = "numeric",
    method = "character", data = "ANY", control = "ControlEL"
  )
)


#' \linkS4class{ELT} class
#'
#' S4 class for empirical likelihood test.
#'
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the (constrained) optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logp A numeric vector of the log probabilities of the (constrained)
#'   empirical likelihood.
#' @slot logl A single numeric of the (constrained) empirical log-likelihood.
#' @slot loglr A single numeric of the (constrained) empirical log-likelihood
#'   ratio.
#' @slot statistic A single numeric of minus twice the (constrained) empirical
#'   log-likelihood ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the chi-square degrees of freedom of the
#'   statistic.
#' @slot pval A single numeric for the (calibrated) \eqn{p}-value of the
#'   statistic.
#' @slot cv A single numeric for the critical value.
#' @slot rhs A numeric vector for the right-hand side of the hypothesis.
#' @slot lhs A numeric matrix for the left-hand side of the hypothesis.
#' @slot alpha A single numeric for the significance level.
#' @slot calibrate A single character for the calibration method used.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases ELT
#' @examples
#' showClass("ELT")
setClass("ELT",
  slots = c(
    optim = "list", logp = "numeric", logl = "numeric", loglr = "numeric",
    statistic = "numeric", df = "integer", pval = "numeric", cv = "numeric",
    rhs = "numeric", lhs = "matrix", alpha = "numeric", calibrate = "character",
    control = "ControlEL"
  )
)


#' \linkS4class{QGLM} class
#'
#' S4 class for generalized linear models with quasi-likelihood methods. It
#'   inherits from \linkS4class{GLM} class.
#'
#' @details The overall test involves a constrained optimization problem. All
#'   the parameters except for the intercept are constrained to zero. The
#'   `optim` slot contains the results. When there is no intercept, all
#'   parameters are set to zero, and the results need to be understood in terms
#'   of \linkS4class{EL} class since no constrained optimization is involved.
#'   Once the solution is found, the log probabilities (`logp`) and the
#'   (constrained) empirical likelihood values (`logl`, `loglr`, `statistic`)
#'   readily follow, along with the degrees of freedom (`df`) and the
#'   \eqn{p}-value (`pval`). The significance tests for each parameter also
#'   involve constrained optimization problems where only one parameter is
#'   constrained to zero. The `sigTests` slot contains the results.
#' @slot family A [`family`] object used.
#' @slot dispersion A single numeric for the estimated dispersion parameter.
#' @slot sigTests A list of the following results of significance tests:
#'   * `statistic` A numeric vector of minus twice the (constrained) empirical
#'   log-likelihood ratios with asymptotic chi-square distributions.
#'   * `iterations` An integer vector for the number of iterations performed for
#'   each parameter.
#'   * `convergence` A logical vector for the convergence status of each
#'   parameter.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot call A matched call.
#' @slot terms A [`terms`] object used.
#' @slot misc A list of various outputs obtained from the model fitting process.
#'   They are used in other generics and methods.
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the (constrained) optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#' @slot logp A numeric vector of the log probabilities of the (constrained)
#'   empirical likelihood.
#' @slot logl A single numeric of the (constrained) empirical log-likelihood.
#' @slot loglr A single numeric of the (constrained) empirical log-likelihood
#'   ratio.
#' @slot statistic A single numeric of minus twice the (constrained) empirical
#'   log-likelihood ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weights A numeric vector of the re-scaled weights used for the model
#'   fitting.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot data A numeric matrix of the data for the model fitting.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases QGLM
#' @examples
#' showClass("QGLM")
setClass("QGLM", contains = "GLM")


#' \linkS4class{SD} class
#'
#' S4 class for standard deviation. It inherits from \linkS4class{EL} class.
#'
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the specified parameters.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logp A numeric vector of the log probabilities of the empirical
#'   likelihood.
#' @slot logl A single numeric of the empirical log-likelihood.
#' @slot loglr A single numeric of the empirical log-likelihood ratio.
#' @slot statistic A single numeric of minus twice the empirical log-likelihood
#'   ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weights A numeric vector of the re-scaled weights used for the model
#'   fitting.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot data A numeric matrix of the data for the model fitting.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases SD
#' @examples
#' showClass("SD")
setClass("SD", contains = "EL")


#' \linkS4class{SummaryEL} class
#'
#' S4 class for a summary of \linkS4class{EL} objects.
#'
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the specified parameters.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logl A single numeric of the empirical log-likelihood.
#' @slot loglr A single numeric of the empirical log-likelihood ratio.
#' @slot statistic A single numeric of minus twice the empirical log-likelihood
#'   ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weighted A single logical for whether the data are weighted or not.
#' @slot coefficients A numeric vector of the maximum empirical likelihood
#'   estimates of the parameters.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases SummaryEL
#' @examples
#' showClass("SummaryEL")
setClass("SummaryEL", slots = c(
  optim = "list", logl = "numeric", loglr = "numeric", statistic = "numeric",
  df = "integer", pval = "numeric", nobs = "integer", npar = "integer",
  weighted = "logical", coefficients = "numeric", method = "character",
  control = "ControlEL"
))


#' \linkS4class{SummaryELMT} class
#'
#' S4 class for a summary of \linkS4class{ELMT} objects.
#'
#' @slot aliased A named logical vector showing if the original coefficients are
#'   aliased.
#' @aliases SummaryELMT
#' @examples
#' showClass("SummaryELMT")
setClass("SummaryELMT", slots = c(
  estimates = "list", statistic = "numeric", df = "integer", pval = "numeric",
  cv = "numeric", rhs = "numeric", lhs = "matrix", alpha = "numeric",
  calibrate = "character"
))


#' \linkS4class{SummaryELT} class
#'
#' S4 class for a summary of \linkS4class{ELT} objects.
#'
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the (constrained) optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logl A single numeric of the (constrained) empirical log-likelihood.
#' @slot loglr A single numeric of the (constrained) empirical log-likelihood
#'   ratio.
#' @slot statistic A single numeric of minus twice the (constrained) empirical
#'   log-likelihood ratio with an asymptotic chi-square distribution.
#' @slot df A single integer for the chi-square degrees of freedom of the
#'   statistic.
#' @slot pval A single numeric for the (calibrated) \eqn{p}-value of the
#'   statistic.
#' @slot cv A single numeric for the critical value.
#' @slot rhs A numeric vector for the right-hand side of the hypothesis.
#' @slot lhs A numeric matrix for the left-hand side of the hypothesis.
#' @slot alpha A single numeric for the significance level.
#' @slot calibrate A single character for the calibration method used.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases SummaryELT
#' @examples
#' showClass("SummaryELT")
setClass("SummaryELT", slots = c(
  optim = "list", logl = "numeric", loglr = "numeric", statistic = "numeric",
  df = "integer", pval = "numeric", cv = "numeric", rhs = "numeric",
  lhs = "matrix", alpha = "numeric", calibrate = "character",
  control = "ControlEL"
))


#' \linkS4class{SummaryLM} class
#'
#' S4 class for a summary of \linkS4class{LM} objects.
#'
#' @slot coefficients A numeric matrix of the results of significance tests.
#' @slot intercept A single logical for whether the given model has an intercept
#'   term or not.
#' @slot na.action Information returned by [`model.frame`] on the special
#'   handling of `NA`s.
#' @slot call A matched call.
#' @slot terms A [`terms`] object used.
#' @slot aliased A named logical vector showing if the original coefficients are
#'   aliased.
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the (constrained) optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logl A single numeric of the empirical log-likelihood.
#' @slot loglr A single numeric of the empirical log-likelihood ratio.
#' @slot statistic A single numeric of minus twice the (constrained) empirical
#'   log-likelihood ratio for the overall test.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weighted A single logical for whether the data are weighted or not.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases SummaryLM
#' @examples
#' showClass("SummaryLM")
setClass("SummaryLM", slots = c(
  coefficients = "matrix", intercept = "logical", na.action = "ANY",
  call = "call", terms = "terms", aliased = "logical", optim = "list",
  logl = "numeric", loglr = "numeric", statistic = "numeric", df = "integer",
  pval = "numeric", nobs = "integer", npar = "integer", weighted = "logical",
  method = "character", control = "ControlEL"
))


#' \linkS4class{SummaryGLM} class
#'
#' S4 class for a summary of \linkS4class{GLM} objects. It inherits from
#'   \linkS4class{SummaryLM} class.
#'
#' @slot family A [`family`] object used.
#' @slot dispersion A single numeric for the estimated dispersion parameter.
#' @slot coefficients A numeric matrix of the results of significance tests.
#' @slot intercept A single logical for whether the given model has an intercept
#'   term or not.
#' @slot na.action Information returned by [`model.frame`] on the special
#'   handling of `NA`s.
#' @slot call A matched call.
#' @slot terms A [`terms`] object used.
#' @slot aliased A named logical vector showing if the original coefficients are
#'   aliased.
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the (constrained) optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logl A single numeric of the empirical log-likelihood.
#' @slot loglr A single numeric of the empirical log-likelihood ratio.
#' @slot statistic A single numeric of minus twice the (constrained) empirical
#'   log-likelihood ratio for the overall test.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weighted A single logical for whether the data are weighted or not.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases SummaryGLM
#' @examples
#' showClass("SummaryGLM")
setClass("SummaryGLM",
  slots = c(family = "family", dispersion = "numeric"), contains = "SummaryLM"
)


#' \linkS4class{SummaryQGLM} class
#'
#' S4 class for a summary of \linkS4class{QGLM} objects. It inherits from
#'   \linkS4class{SummaryGLM} class.
#'
#' @slot family A [`family`] object used.
#' @slot dispersion A single numeric for the estimated dispersion parameter.
#' @slot coefficients A numeric matrix of the results of significance tests.
#' @slot intercept A single logical for whether the given model has an intercept
#'   term or not.
#' @slot na.action Information returned by [`model.frame`] on the special
#'   handling of `NA`s.
#' @slot call A matched call.
#' @slot terms A [`terms`] object used.
#' @slot aliased A named logical vector showing if the original coefficients are
#'   aliased.
#' @slot optim A list of the following optimization results:
#'   * `par` A numeric vector of the solution to the (constrained) optimization
#'   problem.
#'   * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'   problem corresponding to `par`.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#'   * `cstr` A single logical for whether constrained EL optimization is
#'   performed or not.
#' @slot logl A single numeric of the empirical log-likelihood.
#' @slot loglr A single numeric of the empirical log-likelihood ratio.
#' @slot statistic A single numeric of minus twice the (constrained) empirical
#'   log-likelihood ratio for the overall test.
#' @slot df A single integer for the degrees of freedom of the statistic.
#' @slot pval A single numeric for the \eqn{p}-value of the statistic.
#' @slot nobs A single integer for the number of observations.
#' @slot npar A single integer for the number of parameters.
#' @slot weighted A single logical for whether the data are weighted or not.
#' @slot method A single character for the method dispatch in internal
#'   functions.
#' @slot control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @aliases SummaryQGLM
#' @examples
#' showClass("SummaryQGLM")
setClass("SummaryQGLM", contains = "SummaryGLM")
