#' Describe a hypothesis
#'
#' Describe a hypothesis in [print()] for an object of class \linkS4class{ELT}.
#'
#' @param rhs A numeric vector.
#' @param lhs A numeric matrix.
#' @param pnames A character vector.
#' @param digits A single integer for the number of significant digits to be
#'   passed to [format()].
#' @return A character vector for symbolic description of a hypothesis.
#' @noRd
describe_hypothesis <- function(rhs, lhs, pnames, digits) {
  if (is.null(pnames)) {
    pnames <- if (ncol(lhs) == 1L) "par" else paste0("par", seq_len(ncol(lhs)))
  }
  r <- format.default(rhs, digits = digits)
  out <- vector("character", length = nrow(lhs))
  for (i in seq_len(nrow(lhs))) {
    # Indices of the nonzero elements
    idx <- abs(lhs[i, ]) > sqrt(.Machine$double.eps)
    h <- as.character(lhs[i, idx])
    # Append "+" symbols to the positive elements
    h <- ifelse(h > 0, paste0("+", h), h)
    # Add a space after the symbols
    h <- gsub("-", "- ", h, fixed = TRUE)
    h <- gsub("+", "+ ", h, fixed = TRUE)
    # Replace any "+ 1" ("- 1") with "+" ("-")
    h <- gsub("^\\+\\s1$", "+", h)
    h <- gsub("^\\-\\s1$", "-", h)
    # Append the names of the corresponding variables and concatenate the vector
    h <- paste(h, pnames[idx], collapse = " ")
    # Remove any leading " + "
    h <- gsub("^\\+\\s", "", h)
    # Replace any leading "- " with "-"
    h <- gsub("^\\-\\s", "-", h)
    # Append the corresponding element of `rhs`
    h <- paste(h, "=", r[i])
    out[i] <- h
  }
  out
}
