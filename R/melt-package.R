## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib melt, .registration = TRUE
## usethis namespace: end
#' @importFrom graphics points polygon plot.default text
#' @importFrom methods getDataPart is new
#' @importFrom stats .getXlevels fitted gaussian glm.fit is.empty.model lm.fit
#'   lm.wfit model.extract model.matrix model.offset model.response
#'   model.weights naprint pchisq pf printCoefmat qchisq qf quantile setNames
#' @importFrom utils head tail
#' @references
#'   Kim E, MacEachern SN, Peruggia M (2023).
#'   ``Empirical likelihood for the analysis of experimental designs.''
#'   \emph{Journal of Nonparametric Statistics}, **35**(4), 709--732.
#'   \doi{10.1080/10485252.2023.2206919}.
#' @keywords internal
"_PACKAGE"
