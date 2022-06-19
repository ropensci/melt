% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/coef-methods.R
\name{coef}
\alias{coef}
\alias{coef,EL-method}
\title{Model coefficients}
\usage{
\S4method{coef}{EL}(object, ...)
}
\arguments{
\item{object}{A fitted \linkS4class{EL} object.}

\item{...}{Not used.}
}
\value{
A numeric vector of the maximum empirical likelihood estimates.
}
\description{
Extracts maximum empirical likelihood estimates from a model.
}
\examples{
fit <- el_lm(formula = mpg ~ wt, data = mtcars)
coef(fit)
}