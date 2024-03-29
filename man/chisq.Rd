% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/chisq-methods.R
\name{chisq}
\alias{chisq}
\alias{chisq,EL-method}
\alias{chisq,ELMT-method}
\alias{chisq,ELT-method}
\alias{chisq,SummaryEL-method}
\alias{chisq,SummaryELMT-method}
\alias{chisq,SummaryELT-method}
\alias{chisq,SummaryLM-method}
\title{Chi-square statistic}
\usage{
\S4method{chisq}{EL}(object, ...)

\S4method{chisq}{ELMT}(object, ...)

\S4method{chisq}{ELT}(object, ...)

\S4method{chisq}{SummaryEL}(object, ...)

\S4method{chisq}{SummaryELMT}(object, ...)

\S4method{chisq}{SummaryELT}(object, ...)

\S4method{chisq}{SummaryLM}(object, ...)
}
\arguments{
\item{object}{An object that contains the chi-square statistic.}

\item{...}{Further arguments passed to methods.}
}
\value{
The form of the value returned by \code{\link[=chisq]{chisq()}} depends on the class of
its argument.
}
\description{
Extracts the chi-square statistic from a model.
}
\section{Methods (by class)}{
\itemize{
\item \code{chisq(EL)}: Extracts the chi-square statistic.

\item \code{chisq(ELMT)}: Extracts the vector of chi-square statistics.

\item \code{chisq(ELT)}: Extracts the chi-square statistic.

\item \code{chisq(SummaryEL)}: Extracts the chi-square statistic.

\item \code{chisq(SummaryELMT)}: Extracts the vector of chi-square statistics.

\item \code{chisq(SummaryELT)}: Extracts the chi-square statistic.

\item \code{chisq(SummaryLM)}: Extracts the chi-square statistic for the overall test of
the model.

}}
\examples{
data("precip")
fit <- el_mean(precip, par = 40)
chisq(fit)
}
\seealso{
\linkS4class{EL}, \linkS4class{ELMT}, \linkS4class{ELT}, \code{\link[=pVal]{pVal()}}
}
