check_control <- function(control = list()) {
  ctrl <- list(maxit = 100L, abstol = 1e-06, threshold = NULL)
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

  # abstol: Numeric (positive, finite)
  ctrl$abstol <- tryCatch(as.numeric(ctrl$abstol),
                          warning = function(w) NA, error = function(e) NA)
  if (any(length(ctrl$abstol) != 1L, is.na(ctrl$abstol),
          is.infinite(ctrl$abstol)))
    stop("'abstol' is not a number")
  if (ctrl$abstol < .Machine$double.eps)
    stop("'abstol' is not a positive number")

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
