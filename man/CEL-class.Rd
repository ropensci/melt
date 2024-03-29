% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{CEL-class}
\alias{CEL-class}
\alias{CEL}
\title{\linkS4class{CEL} class}
\description{
S4 class for constrained empirical likelihood. It inherits from
\linkS4class{EL} class. Note that the \code{optim} slot has constrained
optimization results with respect to the parameters, not the Lagrange
multiplier.
}
\details{
Let \eqn{l(\theta)} denote minus twice the empirical log-likelihood
ratio function. We consider a linear hypothesis of the form
\deqn{L\theta = r,} where the left-hand-side \eqn{L} is a \eqn{q} by
\eqn{p} matrix and the right-hand-side \eqn{r} is a \eqn{q}-dimensional
vector. Under some regularity conditions, \eqn{l(\theta)} converges in
distribution to \eqn{\chi^2_q} under the constraint of hypothesis, i.e.,
\deqn{\min_{\theta: L\theta = r} l(\theta) \to_d \chi^2_q .}

Minimization of \eqn{l(\theta)} with respect to \eqn{\theta} is
computationally expensive since it implicitly involves the
evaluation step as described in \linkS4class{EL}. Further, depending on the
form of \eqn{g(X_i, \theta)} and the constraint, the optimization problem
can be nonconvex and have multiple local minima. For this reason, the
package \pkg{melt} only considers linear hypotheses and performs local
minimization of \eqn{l(\theta)} using projected gradient descent method.
With the orthogonal projection matrix \eqn{P} and a step size \eqn{\gamma},
the algorithm updates \eqn{\theta} as
\deqn{\theta^{(k + 1)} \leftarrow \theta^{(k)} -
  \gamma P \nabla l(\theta^{(k)}),}
where \eqn{\nabla l(\theta^{(k)})} denotes the gradient of \eqn{l} at
\eqn{\theta^{(k)}}. The first order optimality condition is
\eqn{P \nabla l(\theta) = 0}, which is used as the stopping criterion.
}
\section{Slots}{

\describe{
\item{\code{optim}}{A list of the following optimization results:
\itemize{
\item \code{par} A numeric vector of the solution to the constrained optimization
problem.
\item \code{lambda} A numeric vector of the Lagrange multipliers of the dual
problem corresponding to \code{par}.
\item \code{iterations} A single integer for the number of iterations performed.
\item \code{convergence} A single logical for the convergence status.
\item \code{cstr} A single logical for whether constrained EL optimization is
performed or not.
}}

\item{\code{logp}}{A numeric vector of the log probabilities of the constrained
empirical likelihood.}

\item{\code{logl}}{A single numeric of the constrained empirical log-likelihood.}

\item{\code{loglr}}{A single numeric of the constrained empirical log-likelihood
ratio.}

\item{\code{statistic}}{A single numeric of minus twice the constrained empirical
log-likelihood ratio with an asymptotic chi-square distribution.}

\item{\code{df}}{A single integer for the degrees of freedom of the statistic.}

\item{\code{pval}}{A single numeric for the \eqn{p}-value of the statistic.}

\item{\code{nobs}}{A single integer for the number of observations.}

\item{\code{npar}}{A single integer for the number of parameters.}

\item{\code{weights}}{A numeric vector of the re-scaled weights used for the model
fitting.}

\item{\code{coefficients}}{A numeric vector of the maximum empirical likelihood
estimates of the parameters.}

\item{\code{method}}{A single character for the method dispatch in internal
functions.}

\item{\code{data}}{A numeric matrix of the data for the model fitting.}

\item{\code{control}}{An object of class \linkS4class{ControlEL} constructed by
\code{\link[=el_control]{el_control()}}.}
}}

\examples{
showClass("CEL")
}
\references{
Adimari G, Guolo A (2010).
``A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.''
\emph{Statistical Methods & Applications}, \strong{19}(4), 463--476.
\doi{10.1007/s10260-010-0137-9}.

Qin J, Lawless J (1995).
``Estimating Equations, Empirical Likelihood and Constraints on
Parameters.'' \emph{Canadian Journal of Statistics}, \strong{23}(2), 145--159.
\doi{10.2307/3315441}.
}
