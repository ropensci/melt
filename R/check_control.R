check_control <- function(control = list()) {
  ctrl <- list(maxit = 100L, tol = 1e-06, th = NULL)
  nctrl <- names(ctrl)
  ctrl[(ncontrol <- names(control))] <- control
  if (length(nomatch <- ncontrol[!ncontrol %in% nctrl]))
    warning("unknown names in control: ", paste(nomatch, collapse = ", "))

  ##
  # deprecated (threshold, abstol)
  if ("threshold" %in% names(control)) {
    warning("'threshold' is deprecated; please use 'th' instead.",
            call. = FALSE)
  }
  if ("abstol" %in% names(control)) {
    warning("'abstol' is deprecated; please use 'tol' instead.",
            call. = FALSE)
  }
  ##

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

  # th: Numeric (positive, finite)
  if (!is.null(ctrl$th)) {
    ctrl$th <- tryCatch(as.numeric(ctrl$th),
                               warning = function(w) NA, error = function(e) NA)
    if (any(length(ctrl$th) != 1L, is.na(ctrl$th),
            is.infinite(ctrl$th)))
      stop("'th' is not a number")
    if (ctrl$th < .Machine$double.eps)
      stop("'th' is not a positive number")
  }
  ctrl
}
