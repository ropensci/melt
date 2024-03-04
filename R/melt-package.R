## usethis namespace: start
#' @importFrom checkmate assert_class
#' @importFrom checkmate assert_int
#' @importFrom checkmate assert_logical
#' @importFrom checkmate assert_matrix
#' @importFrom checkmate assert_number
#' @importFrom checkmate assert_numeric
#' @importFrom checkmate assert_string
#' @importFrom checkmate test_numeric
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
#' @references
#'   Kim E, MacEachern SN, Peruggia M (2024).
#'   ``melt: Multiple Empirical Likelihood Tests in R.''
#'   \emph{Journal of Statistical Software}, **108**(5), 1--33.
#'   \doi{10.18637/jss.v108.i05}.
#' @keywords internal
"_PACKAGE"
