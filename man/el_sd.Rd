% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/el_sd.R
\name{el_sd}
\alias{el_sd}
\title{Empirical likelihood for the standard deviation}
\usage{
el_sd(x, mean, sd, weights = NULL, control = el_control())
}
\arguments{
\item{x}{A numeric vector, or an object that can be coerced to a numeric
vector.}

\item{mean}{A single numeric for the (known) mean value.}

\item{sd}{A positive single numeric for the parameter value to be tested.}

\item{weights}{An optional numeric vector of weights to be used in the
fitting process. The length of the vector must be the same as the length of
\code{x}. Defaults to \code{NULL}, corresponding to identical weights. If non-\code{NULL},
weighted empirical likelihood is computed.}

\item{control}{An object of class \linkS4class{ControlEL} constructed by
\code{\link[=el_control]{el_control()}}.}
}
\value{
An object of class \linkS4class{SD}.
}
\description{
Computes empirical likelihood for the standard deviation.
}
\details{
Let \eqn{X_i} be independent and identically random variable from an
unknown distribution \eqn{P} for \eqn{i = 1, \dots, n}. We assume that
\eqn{{\textrm{E}[X_i]} = {\mu_0}} is known and that \eqn{P} has a variance
\eqn{\sigma_0^2}. Given a value of \eqn{\sigma}, the
(profile) empirical likelihood ratio is defined by
\deqn{R(\sigma) =
  \max_{p_i}\left\{\prod_{i = 1}^n np_i :
  \sum_{i = 1}^n p_i (X_i - \mu_0)^2 = \sigma^2,\
  p_i \geq 0,\
  \sum_{i = 1}^n p_i = 1
  \right\}.}
\code{\link[=el_sd]{el_sd()}} computes the empirical log-likelihood ratio statistic
\eqn{-2\log R(\sigma)}, along with other values in \linkS4class{SD}.
}
\examples{
data("women")
x <- women$height
w <- women$weight
fit <- el_sd(x, mean = 65, sd = 5, weights = w)
fit
summary(fit)
}
\seealso{
\linkS4class{EL}, \linkS4class{SD}, \code{\link[=el_mean]{el_mean()}}, \code{\link[=elt]{elt()}},
\code{\link[=el_control]{el_control()}}
}
