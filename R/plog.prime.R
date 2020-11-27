#' @export
plog.prime <- function(x, threshold) {
  if(anyNA(x) == T) stop("NA input found in plog_prime function")
  output <- vector(mode = "numeric", length = length(x))
  test <- drop(x >= threshold)
  output[test] <- x[test]^-1
  output[!test] <- 2 * threshold^-1 - x[!test] / threshold^2
  output
}
