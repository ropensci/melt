#' Chi-square statistic
#'
#' Extracts the chi-square statistic from a model.
#'
#' @param object An object that contains the chi-square statistic.
#' @param ... Further arguments passed to methods.
#' @return The form of the value returned by [chisq()] depends on the class of
#'   its argument.
#' @seealso \linkS4class{EL}, \linkS4class{ELMT}, \linkS4class{ELT}, [pVal()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' chisq(fit)
#' @exportMethod chisq
setGeneric("chisq", function(object, ...) standardGeneric("chisq"))


#' Model coefficients
#'
#' Extracts the maximum empirical likelihood estimates from a model.
#'
#' @param object An object that contains the maximum empirical likelihood
#'   estimates.
#' @param ... Further arguments passed to methods.
#' @return The form of the value returned by [coef()] depends on the class of
#'   its argument.
#' @seealso \linkS4class{EL}, \linkS4class{ELMT}
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ wt, data = mtcars)
#' coef(fit)
#' @exportMethod coef
setGeneric("coef", function(object, ...) standardGeneric("coef"))


#' Confidence interval for model parameters
#'
#' Computes confidence intervals for one or more parameters in a model.
#'
#' @param object An object that inherits from \linkS4class{EL} or
#'   \linkS4class{ELMT}.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing, all
#'   parameters are considered.
#' @param level A single numeric for the confidence level required. Defaults to
#'   `0.95`.
#' @param cv A single numeric for the critical value for calibration of
#'   empirical likelihood ratio statistic. Defaults to `NULL` and set to
#'   `qchisq(level, 1L)`. If non-`NULL`, `level` is ignored.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()]. Defaults to `NULL` and inherits the `control` slot in
#'   `object`.
#' @return A matrix with columns giving lower and upper confidence limits for
#'  each parameter. In contrast to other methods that rely on studentization,
#'  the lower and upper limits obtained from empirical likelihood do not
#'  correspond to the `(1 - level) / 2` and `1 - (1 - level) / 2` in %,
#'  respectively.
#' @references Owen A (1990).
#'   “Empirical Likelihood Ratio Confidence Regions.”
#'   \emph{The Annals of Statistics}, 18(1), 90--120.
#'   \doi{10.1214/aos/1176347494}.
#' @seealso \linkS4class{EL}, \linkS4class{ELMT}, [confreg()], [elt()],
#'   [el_control()]
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' confint(fit, parm = c(2, 3))
#' @exportMethod confint
setGeneric("confint", function(object, parm, level = 0.95, ...)
  standardGeneric("confint")
)


#' Confidence region for model parameters
#'
#' Computes boundary points of a two-dimensional confidence region for model
#'   parameters.
#'
#' @param object An object that inherits from \linkS4class{EL}.
#' @param parm A specification of which parameters are to be given a confidence
#'   region, either a vector of numbers or a vector of names. It must be a
#'   vector of length two of the form `c(x, y)`. If missing, the first two
#'   parameter in `object` are considered.
#' @param level A single numeric for the confidence level required. Defaults to
#'   `0.95`. It is ignored if `cv` is non-`NULL`.
#' @param cv A single numeric for the critical value for calibration of
#'   empirical likelihood ratio statistic. Defaults to NULL and set to
#'   `qchisq(level, 2L)`. It must be compatible with the `th` value in
#'   `control`.
#' @param npoints A single integer for the number of boundary points to compute.
#'   Defaults to `50`.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()]. Defaults to `NULL` and inherits the `control` slot in
#'   `object`.
#' @return An object of class \linkS4class{ConfregEL}.
#' @references Owen A (1990).
#'   “Empirical Likelihood Ratio Confidence Regions.”
#'   \emph{The Annals of Statistics}, 18(1), 90--120.
#'   \doi{10.1214/aos/1176347494}.
#' @seealso \linkS4class{EL}, [confint()], [elt()], [plot()], [el_control()]
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ wt + qsec, data = mtcars)
#' cr <- confreg(fit, parm = c(2, 3), cv = qchisq(0.90, 2))
#' plot(cr)
#' @exportMethod confreg
setGeneric("confreg", function(object,
                               parm,
                               level = 0.95,
                               cv = NULL,
                               npoints = 50L,
                               control = NULL) {
  standardGeneric("confreg")
})


#' Convergence check
#'
#' Extracts the convergence status from a model.
#'
#' @param object An object that contains the convergence status.
#' @param ... Further arguments passed to methods.
#' @return A single logical.
#' @seealso \linkS4class{CEL}, \linkS4class{EL}, \linkS4class{ELT}, [getOptim()]
#' @usage NULL
#' @examples
#' ## Convergence check for the overall model test
#' data("mtcars")
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' conv(fit)
#' @exportMethod conv
setGeneric("conv", function(object, ...) standardGeneric("conv"))


#' Critical value
#'
#' Extracts the critical value from a model.
#'
#' @param object An object that contains the critical value.
#' @param ... Further arguments passed to methods.
#' @return A single numeric.
#' @seealso \linkS4class{ELMT}, \linkS4class{ELT}
#' @usage NULL
#' @examples
#' ## F-calibrated critical value
#' data("precip")
#' fit <- el_mean(precip, 30)
#' elt <- elt(fit, rhs = 34, calibrate = "f")
#' critVal(elt)
#' @exportMethod critVal
setGeneric("critVal", function(object, ...) standardGeneric("critVal"))


#' Empirical likelihood displacement
#'
#' Computes empirical likelihood displacement for model diagnostics and outlier
#'   detection.
#'
#' @param object An object that inherits from \linkS4class{EL}.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()]. Defaults to `NULL` and inherits the `control` slot in
#'   `object`.
#' @details Let \eqn{L(\theta)} be the empirical log-likelihood function based
#'   on the full sample with \eqn{n} observations. The maximum empirical
#'   likelihood estimate is denoted by \eqn{\hat{\theta}}. Consider a reduced
#'   sample with the \eqn{i}th observation deleted and the corresponding
#'   estimate \eqn{\hat{\theta}_{(i)}}. The empirical likelihood displacement is
#'   defined by
#'   \deqn{\textrm{ELD}_i = 2\{L(\hat{\theta}) - L(\hat{\theta}_{(i)})\}.}
#'   If \eqn{\textrm{ELD}_i } is large, then the \eqn{i}th observation is an
#'   influential point and can be inspected as a possible outlier. `eld`
#'   computes \eqn{\textrm{ELD}_i } for \eqn{i = 1, \dots, n }.
#' @return An object of class \linkS4class{ELD}.
#' @references Lazar NA (2005).
#'   “Assessing the Effect of Individual Data Points on Inference From Empirical
#'   Likelihood.”
#'   \emph{Journal of Computational and Graphical Statistics}, 14(3), 626–642.
#'   \doi{10.1198/106186005X59568}.
#' @references Zhu H, Ibrahim JG, Tang N, Zhang H (2008).
#'   “Diagnostic Measures for Empirical Likelihood of General Estimating
#'   Equations.” \emph{Biometrika}, 95(2), 489--507.
#'   \doi{10.1093/biomet/asm094}.
#' @seealso \linkS4class{EL}, \linkS4class{ELD}, [el_control()], [plot()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 30)
#' eld <- eld(fit)
#' plot(eld)
#' @exportMethod eld
setGeneric("eld", function(object, control = NULL) {
  standardGeneric("eld")
})


#' Empirical likelihood multiple tests
#'
#' Tests multiple linear hypotheses simultaneously.
#'
#' @param object An object that inherits from \linkS4class{EL}.
#' @param rhs A numeric vector (column matrix) or a list of numeric vectors for
#'   the right-hand sides of hypotheses. Defaults to `NULL`. See ‘Details’.
#' @param lhs A list or a numeric matrix for the left-hand sides of hypotheses.
#'   For a list `lhs`, each element must be specified as a single instance of
#'   the `lhs` in [elt()]. For a matrix `lhs`, each row gives a linear
#'   combination of the parameters in `object`. The number of columns must be
#'   equal to the number of parameters. Defaults to `NULL`. See ‘Details’.
#' @param alpha A single numeric for the overall significance level. Defaults to
#'   `0.05`.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()]. Defaults to `NULL` and inherits the `control` slot in
#'   `object`.
#' @details [elmt()] tests multiple hypotheses simultaneously. Each hypothesis
#'   corresponds to the constrained empirical likelihood ratio described in
#'   \linkS4class{CEL}. `rhs` and `lhs` cannot be both `NULL`. The right-hand
#'   side and left-hand side of each hypothesis must be specified as described
#'   in [elt()].
#'
#'   For specifying linear contrasts more conveniently, `rhs` and `lhs` also
#'   take a numeric vector and a numeric matrix, respectively. Each element of
#'   `rhs` and each row of `lhs` correspond to a contrast (hypothesis).
#'
#'   The vector of empirical likelihood ratio statistics asymptotically follows
#'   a multivariate chi-square distribution under the complete null hypothesis.
#'   The multiple testing procedure asymptotically controls the family-wise
#'   error rate at the level `alpha`. Based on the distribution of the maximum
#'   of the test statistics, the adjusted p-values are estimated by Monte Carlo
#'   simulation.
#' @return An object of class of \linkS4class{ELMT}.
#' @references Kim E, MacEachern S, Peruggia M (2021).
#'   “Empirical Likelihood for the Analysis of Experimental Designs.”
#'   arxiv:2112.09206. URL <https://arxiv.org/abs/2112.09206>.
#' @seealso \linkS4class{EL}, \linkS4class{ELMT}, [elt()], [el_control()]
#' @usage NULL
#' @examples
#' ## Bivariate mean (list `rhs` & no `lhs`)
#' set.seed(143)
#' data("women")
#' fit <- el_mean(women, par = c(65, 135))
#' rhs <- list(c(64, 133), c(66, 140))
#' elmt(fit, rhs = rhs)
#'
#' ## Pairwise comparison (no `rhs` & list `lhs`)
#' data("clothianidin")
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' lhs2 <- list(
#'   "trtNaked - trtFungicide",
#'   "trtFungicide - trtLow",
#'   "trtLow - trtHigh"
#' )
#' elmt(fit2, lhs = lhs2)
#'
#' ## Arbitrary hypotheses (list `rhs` & list `lhs`)
#' data("mtcars")
#' fit3 <- el_lm(mpg ~ wt + qsec, data = mtcars)
#' lhs3 <- list(c(1, 4, 0), rbind(c(0, 1, 0), c(0, 0, 1)))
#' rhs3 <- list(0, c(-6, 1))
#' elmt(fit3, rhs = rhs3, lhs = lhs3)
#' @exportMethod elmt
setGeneric("elmt", function(object,
                            rhs = NULL,
                            lhs = NULL,
                            alpha = 0.05,
                            control = NULL) {
  standardGeneric("elmt")
})


#' Empirical likelihood test
#'
#' Tests a linear hypothesis.
#'
#' @param object An object that inherits from \linkS4class{EL}.
#' @param rhs A numeric vector or a column matrix for the right-hand side of
#'   hypothesis, with as many entries as the rows in `lhs`. Defaults to `NULL`.
#'   See ‘Details’.
#' @param lhs A numeric matrix or a vector (treated as a row matrix) for the
#'   left-hand side of a hypothesis. Each row gives a linear combination of the
#'   parameters in `object`. The number of columns must be equal to the number
#'   of parameters. Or a character vector with a symbolic description of the
#'   hypothesis is allowed. Defaults to `NULL`. See ‘Details’.
#' @param alpha A single numeric for the significance level. Defaults to `0.05`.
#' @param calibrate A single character for the calibration method. It is
#'   case-insensitive and must be one of `"chisq"`, `"boot"`, or `"f"`.
#'   Defaults to `"chisq"`. See ‘Details’.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()]. Defaults to `NULL` and inherits the `control` slot in
#'   `object`.
#' @details [elt()] performs the constrained minimization of \eqn{l(\theta)}
#'   described in \linkS4class{CEL}. `rhs` and `lhs` cannot be both `NULL`. For
#'   non-`NULL` `lhs`, it is required that `lhs` have full row rank
#'   \eqn{q \leq p} and \eqn{p} be equal to the number of parameters in the
#'   `object`.
#'
#'   Depending on the specification of `rhs` and `lhs`, we have the following
#'   three cases:
#'   1. If both `rhs` and `lhs` are non-`NULL`, the constrained minimization
#'   is performed with the right-hand side \eqn{r} and the left-hand side
#'   \eqn{L} as
#'   \deqn{\inf_{\theta: L\theta = r} l(\theta).}
#'   1. If `rhs` is `NULL`, \eqn{r} is set to the zero vector as
#'   \eqn{\inf_{\theta: L\theta = 0} l(\theta)}.
#'   1. If `lhs` is `NULL`, \eqn{L} is set to the identity matrix and the
#'   problem reduces to evaluating at \eqn{r} as \eqn{l(r)}.
#'
#'   `calibrate` specifies the calibration method used. Three methods are
#'   available: `"chisq"` (chi-square calibration), `"boot"` (bootstrap
#'   calibration), and `"f"` (\eqn{F} calibration). `"boot"` is applicable only
#'   when `lhs` is `NULL`. The `nthreads`, `seed`, and `B` slots in `control`
#'   apply to the bootstrap procedure. `"f"` is applicable only to the mean
#'   parameter when `lhs` is `NULL`.
#' @return An object of class of \linkS4class{ELT}. If `lhs` is non-`NULL`, the
#'   `optim` slot corresponds to that of \linkS4class{CEL}. Otherwise, it
#'   corresponds to that of \linkS4class{EL}.
#' @references Adimari G, Guolo A (2010).
#'   “A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.”
#'   \emph{Statistical Methods & Applications}, 19(4), 463--476.
#'   \doi{10.1007/s10260-010-0137-9}.
#' @references Qin J, Lawless J (1995).
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   \emph{Canadian Journal of Statistics}, 23(2), 145--159.
#'   \doi{10.2307/3315441}.
#' @seealso \linkS4class{EL}, \linkS4class{ELT}, [elmt()], [el_control()]
#' @usage NULL
#' @examples
#' ## F calibration for the mean
#' data("precip")
#' fit <- el_mean(precip, 32)
#' elt(fit, rhs = 32, calibrate = "f")
#'
#' ## Test of no treatment effect
#' data("clothianidin")
#' contrast <- matrix(c(
#'   1, -1, 0, 0,
#'   0, 1, -1, 0,
#'   0, 0, 1, -1
#' ), byrow = TRUE, nrow = 3)
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' elt(fit2, lhs = contrast)
#'
#' ## A symbolic description of the same hypothesis
#' elt(fit2, lhs = c(
#'   "trtNaked - trtFungicide",
#'   "trtFungicide - trtLow",
#'   "trtLow - trtHigh"
#' ))
#' @exportMethod elt
setGeneric("elt", function(object,
                           rhs = NULL,
                           lhs = NULL,
                           alpha = 0.05,
                           calibrate = "chisq",
                           control = NULL) {
  standardGeneric("elt")
})


#' Degrees of freedom
#'
#' Extracts the degrees of freedom from a model.
#'
#' @param object An object that contains the degrees of freedom.
#' @return An integer vector.
#' @seealso \linkS4class{EL}, \linkS4class{ELMT}, \linkS4class{ELT}
#' @usage NULL
#' @examples
#' data("faithful")
#' fit <- el_mean(faithful, par = c(3.5, 70))
#' getDF(fit)
#' @exportMethod getDF
setGeneric("getDF", function(object) standardGeneric("getDF"))


#' Optimization results
#'
#' Extracts the optimization results from a model.
#'
#' @param object An object that contains the optimization results.
#' @param ... Further arguments passed to methods.
#' @return A list with the following optimization results:
#'   * `par` A numeric vector of the parameter value. See the documentation of
#'   \linkS4class{EL} and \linkS4class{CEL}.
#'   * `lambda` A numeric vector of the Lagrange multipliers.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#' @seealso \linkS4class{EL}, \linkS4class{ELT}, [sigTests()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' getOptim(fit)
#' @exportMethod getOptim
setGeneric("getOptim", function(object, ...) standardGeneric("getOptim"))


#' Empirical log-likelihood
#'
#' Extracts the empirical log-likelihood from a model.
#'
#' @param object An object that contains the empirical log-likelihood.
#' @param ... Further arguments passed to methods.
#' @return A single numeric.
#' @references Baggerly KA (1998).
#'   “Empirical Likelihood as a Goodness-of-Fit Measure.” \emph{Biometrika},
#'   85(3), 535--547. \doi{10.1093/biomet/asm094}.
#' @seealso \linkS4class{EL}, \linkS4class{ELT}
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' logL(fit)
#' @exportMethod logL
setGeneric("logL", function(object, ...) standardGeneric("logL"))


#' Maximum empirical log-likelihood
#'
#' @description Extracts empirical log-likelihood from a model evaluated at the
#'   estimated coefficients.
#'
#' \emph{This function is deprecated and will be removed in a future release.}
#'
#' @param object An object that inherits from \linkS4class{EL}.
#' @param ... Further arguments passed to methods.
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
#'   The maximum empirical likelihood estimator \eqn{\hat{\theta}} solves
#'   \eqn{n^{-1}\sum_{i = 1}^n g(X_i, \hat{\theta}) = 0} and yields
#'   \eqn{p_i = 1/n} for \eqn{i = 1, \dots, n}. [logLik()] gives \eqn{-n\log n},
#'   the maximum empirical log-likelihood. Use [logL()] instead to extract the
#'   (constrained) empirical log-likelihood computed from a model.
#' @return An object of class \linkS4class{logLikEL}.
#' @seealso \linkS4class{EL}, [logL()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' logLik(fit)
#' @exportMethod logLik
setGeneric("logLik", function(object, ...) standardGeneric("logLik"))


#' Empirical log-likelihood ratio
#'
#' Extracts the empirical log-likelihood ratio from a model.
#'
#' @param object An object that contains the empirical log-likelihood ratio.
#' @param ... Further arguments passed to methods.
#' @return A single numeric.
#' @references Baggerly KA (1998).
#'   “Empirical Likelihood as a Goodness-of-Fit Measure.” \emph{Biometrika},
#'   85(3), 535--547. \doi{10.1093/biomet/asm094}.
#' @seealso \linkS4class{EL}, \linkS4class{ELT}
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' logLR(fit)
#' @exportMethod logLR
setGeneric("logLR", function(object, ...) standardGeneric("logLR"))


#' Log probabilities
#'
#' Extracts log probabilities of empirical likelihood from a model.
#'
#' @param object An object that inherits from \linkS4class{EL} or
#'   \linkS4class{ELT}.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector.
#' @seealso \linkS4class{EL}, \linkS4class{ELT}
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' logProb(fit)
#' @exportMethod logProb
setGeneric("logProb", function(object, ...) standardGeneric("logProb"))


#' Number of observations in a model
#'
#' Extracts the number of observations from a model.
#'
#' @param object An object that contains the number of observations.
#' @param ... Further arguments passed to methods.
#' @return A single integer.
#' @seealso \linkS4class{EL}
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' nobs(fit)
#' @exportMethod nobs
setGeneric("nobs", function(object, ...) standardGeneric("nobs"))


#' Plot methods
#'
#' Provides plot methods for objects.
#'
#' @param x An object to be plotted.
#' @param y Not used.
#' @param ... Further graphical parameters (see [`par`]).
#' @return No return value, called for side effects.
#' @seealso \linkS4class{ConfregEL}, \linkS4class{EL}, \linkS4class{ELD},
#'   [confreg()], [eld()]
#' @usage NULL
#' @examples
#' ## Model
#' data("mtcars")
#' fit <- el_lm(hp ~ wt, data = mtcars)
#'
#' ## Confidence region
#' out1 <- confreg(fit, npoints = 500)
#' plot(out1)
#'
#' ## Empirical likelihood displacement
#' out2 <- eld(fit)
#' plot(out2)
#'
#' ## A shortcut to `ELD`
#' plot(fit)
#' @exportMethod plot
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))


#' Print methods
#'
#' Provides print methods for objects.
#'
#' @param x An object to be printed.
#' @param digits A single integer for the number of significant digits to be
#'   passed to [format()].
#' @param signif.stars A single logical. If `TRUE`, ‘significance stars’ are
#'   printed for each parameter.
#' @param ... Further arguments passed to methods.
#' @return The argument `x` (invisibly).
#' @seealso \linkS4class{EL}, \linkS4class{ELMT}, \linkS4class{ELT},
#'   \linkS4class{LM}
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' print(fit)
#' @exportMethod print
setGeneric("print", function(x, ...) standardGeneric("print"))


#' \eqn{p}-value
#'
#' Extracts the \eqn{p}-value from a model.
#'
#' @param object An object that contains the \eqn{p}-value.
#' @param ... Further arguments passed to methods.
#' @return The form of the value returned by [pVal()] depends on the class of
#'   its argument.
#' @seealso \linkS4class{EL}, \linkS4class{ELMT}, \linkS4class{ELT}, [chisq()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' pVal(fit)
#' @exportMethod pVal
setGeneric("pVal", function(object, ...) standardGeneric("pVal"))


#' Significance tests
#'
#' Extracts the results of significance tests from a model.
#'
#' @param object An object that inherits from \linkS4class{LM} or
#'   \linkS4class{SummaryLM}.
#' @param ... Further arguments passed to methods.
#' @return The form of the value returned by [sigTests()] depends on the
#'   class of its argument.
#' @seealso \linkS4class{LM}, \linkS4class{SummaryLM}, [getOptim()]
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' sigTests(fit)
#' sigTests(summary(fit))
#' @exportMethod sigTests
setGeneric("sigTests", function(object, ...) standardGeneric("sigTests"))


#' Summary methods
#'
#' Provides summary methods for objects.
#'
#' @param object An object for which a summary is desired.
#' @param ... Further arguments passed to methods.
#' @return The form of the value returned by [summary()] depends on the class of
#'   its argument.
#' @seealso \linkS4class{EL}, \linkS4class{ELMT}, \linkS4class{ELT},
#'   \linkS4class{GLM}, \linkS4class{LM}, \linkS4class{QGLM},
#' @usage NULL
#' @examples
#' data("faithful")
#' fit <- el_mean(faithful, par = c(3.5, 70))
#' summary(fit)
#'
#' data("mtcars")
#' fit2 <- el_lm(mpg ~ wt, data = mtcars)
#' summary(fit2)
#' @exportMethod summary
setGeneric("summary", function(object, ...) standardGeneric("summary"))


#' Model weights
#'
#' Extracts weights from model objects. The weights are re-scaled to up to the
#'   total number of observations in the fitting procedure.
#'
#' @param object An object that inherits from \linkS4class{EL}.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector of the re-scaled weights.
#' @references Glenn N, Zhao Y (2007).
#'   “Weighted Empirical Likelihood Estimates and Their Robustness Properties.”
#'   \emph{Computational Statistics & Data Analysis}, 51(10), 5130--5141.
#'   \doi{10.1016/j.csda.2006.07.032}.
#' @seealso \linkS4class{EL}
#' @usage NULL
#' @examples
#' data("airquality")
#' x <- airquality$Wind
#' w <- airquality$Day
#' fit <- el_mean(x, par = 10, weights = w)
#' weights(fit)
#' @exportMethod weights
setGeneric("weights", function(object, ...) standardGeneric("weights"))


setGeneric("getData", function(x) standardGeneric("getData"))
setMethod("getData", "EL", function(x) {
  x@data
})
setMethod("getData", "ELMT", function(x) {
  x@data
})


setGeneric("getEstimates", function(x) standardGeneric("getEstimates"))
setMethod("getEstimates", "EL", function(x) {
  x@coefficients
})
setMethod("getEstimates", "ELMT", function(x) {
  x@estimates
})
setMethod("getEstimates", "QGLM", function(x) {
  c(x@coefficients, setNames(x@dispersion, "phi"))
})
setMethod("getEstimates", "SummaryELMT", function(x) {
  x@estimates
})


setGeneric("getMethodEL", function(x) standardGeneric("getMethodEL"))
setMethod("getMethodEL", "EL", function(x) {
  x@method
})
setMethod("getMethodEL", "ELMT", function(x) {
  x@method
})
setMethod("getMethodEL", "SummaryEL", function(x) {
  x@method
})
setMethod("getMethodEL", "SummaryLM", function(x) {
  x@method
})


setGeneric("getNumPar", function(x) standardGeneric("getNumPar"))
setMethod("getNumPar", "EL", function(x) {
  x@npar
})
setMethod("getNumPar", "SummaryEL", function(x) {
  x@npar
})
setMethod("getNumPar", "SummaryLM", function(x) {
  x@npar
})


setGeneric("getWeights", function(x) standardGeneric("getWeights"))
setMethod("getWeights", "EL", function(x) {
  x@weights
})
setMethod("getWeights", "ELMT", function(x) {
  x@weights
})


setGeneric("getControlEL", function(x) standardGeneric("getControlEL"))
setMethod("getControlEL", "EL", function(x) {
  x@control
})
setMethod("getControlEL", "ELMT", function(x) {
  x@control
})
setMethod("getControlEL", "ELT", function(x) {
  x@control
})
setMethod("getControlEL", "SummaryELT", function(x) {
  x@control
})
setMethod("getControlEL", "SummaryLM", function(x) {
  x@control
})
