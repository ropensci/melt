% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/critVal-methods.R
\name{critVal}
\alias{critVal}
\alias{critVal,ELMT-method}
\alias{critVal,ELT-method}
\alias{critVal,SummaryELMT-method}
\alias{critVal,SummaryELT-method}
\title{Critical value}
\usage{
\S4method{critVal}{ELMT}(object, ...)

\S4method{critVal}{ELT}(object, ...)

\S4method{critVal}{SummaryELMT}(object, ...)

\S4method{critVal}{SummaryELT}(object, ...)
}
\arguments{
\item{object}{An object that contains the critical value.}

\item{...}{Further arguments passed to methods.}
}
\value{
A single numeric.
}
\description{
Extracts the critical value from a model.
}
\examples{
## F-calibrated critical value
data("precip")
fit <- el_mean(precip, 30)
elt <- elt(fit, rhs = 34, calibrate = "f")
critVal(elt)
}
\seealso{
\linkS4class{ELMT}, \linkS4class{ELT}
}
