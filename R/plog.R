#' @noRd
plog <- function(x, threshold) {
  if(anyNA(x) == T) stop("NA input found in plog function")
  output <- vector(mode = "numeric", length = length(x))
  test <- drop(x >= threshold)
  output[test] <- log(x[test])
  output[!test] <- log(threshold) - 1.5 +
    2 * threshold^(-1) * x[!test] - threshold^(-2) * x[!test]^2 / 2
  output
}
