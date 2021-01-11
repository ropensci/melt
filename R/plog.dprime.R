#' @noRd
plog.dprime <- function(x, threshold) {
  if (anyNA(x) == T) stop("NA input found in plog_dprime function")
  output <- vector(mode = "numeric", length = length(x))
  test <- drop(x >= threshold)
  output[test] <- -x[test]^(-2)
  output[!test] <- -(length(x)^2)
  output
}
