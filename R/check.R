check_control <- function(control = list()) {
  ctrl <- list(maxit = 100L, tol = 1e-06, threshold = NULL)
  nctrl <- names(ctrl)
  ctrl[(ncontrol <- names(control))] <- control
  if (length(nomatch <- ncontrol[!ncontrol %in% nctrl]))
    warning("unknown names in control: ", paste(nomatch, collapse = ", "))

  # maxit: Integer (positive)
  ctrl$maxit <- tryCatch(as.integer(ctrl$maxit),
                         warning = function(w) NA, error = function(e) NA)
  if (any(length(ctrl$maxit) != 1L, is.na(ctrl$maxit)))
    stop("'maxit' is not an integer")
  if (ctrl$maxit < 1)
    stop("'maxit' is not a positive integer")

  # tol: Numeric (positive, finite)
  ctrl$tol <- tryCatch(as.numeric(ctrl$tol),
                          warning = function(w) NA, error = function(e) NA)
  if (any(length(ctrl$tol) != 1L, is.na(ctrl$tol),
          is.infinite(ctrl$tol)))
    stop("'tol' is not a number")
  if (ctrl$tol < .Machine$double.eps)
    stop("'tol' is not a positive number")

  # threshold: Numeric (positive, finite)
  if (!is.null(ctrl$threshold)) {
    ctrl$threshold <- tryCatch(as.numeric(ctrl$threshold),
                               warning = function(w) NA, error = function(e) NA)
    if (any(length(ctrl$threshold) != 1L, is.na(ctrl$threshold),
            is.infinite(ctrl$threshold)))
      stop("'threshold' is not a number")
    if (ctrl$threshold < .Machine$double.eps)
      stop("'threshold' is not a positive number")
  }
  ctrl
}
