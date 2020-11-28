#' @noRd
dfp1d <- function(gradient) {
  Q <- diag(1, nrow = 2, ncol = 2)
  A <- c(1, -1)
  LHS <- rbind(cbind(Q, A), cbind(t(A), 0))
  RHS <- -c(gradient, 0)
  unname(solve(LHS, RHS)[1:2])[1]
}
