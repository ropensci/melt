% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/confreg-methods.R
\name{confreg}
\alias{confreg}
\alias{confreg,EL-method}
\title{Confidence region for model parameters}
\usage{
\S4method{confreg}{EL}(object, parm, level = 0.95, cv = NULL, npoints = 50L, control = NULL)
}
\arguments{
\item{object}{An object that inherits from \linkS4class{EL}.}

\item{parm}{A specification of which parameters are to be given a confidence
region, either a vector of numbers or a vector of names. It must be a
vector of length two of the form \code{c(x, y)}. If missing, the first two
parameter in \code{object} are considered.}

\item{level}{A single numeric for the confidence level required. Defaults to
\code{0.95}. It is ignored if \code{cv} is non-\code{NULL}.}

\item{cv}{A single numeric for the critical value for calibration of
empirical likelihood ratio statistic. Defaults to NULL and set to
\code{qchisq(level, 2L)}. It must be compatible with the \code{th} value in
\code{control}.}

\item{npoints}{A single integer for the number of boundary points to compute.
Defaults to \code{50}.}

\item{control}{An object of class \linkS4class{ControlEL} constructed by
\code{\link[=el_control]{el_control()}}. Defaults to \code{NULL} and inherits the \code{control} slot in
\code{object}.}
}
\value{
An object of class \linkS4class{ConfregEL}.
}
\description{
Computes boundary points of a two-dimensional confidence region for model
parameters.
}
\examples{
data("mtcars")
fit <- el_lm(mpg ~ wt + qsec, data = mtcars)
cr <- confreg(fit, parm = c(2, 3), cv = qchisq(0.90, 2))
plot(cr)
}
\references{
Owen A (1990).
``Empirical Likelihood Ratio Confidence Regions.''
\emph{The Annals of Statistics}, \strong{18}(1), 90--120.
\doi{10.1214/aos/1176347494}.
}
\seealso{
\linkS4class{EL}, \code{\link[=confint]{confint()}}, \code{\link[=elt]{elt()}}, \code{\link[=plot]{plot()}}, \code{\link[=el_control]{el_control()}}
}
