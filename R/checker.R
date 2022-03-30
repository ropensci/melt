#' @noRd
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

  # maxit: integer (positive)
  ctrl$maxit <- tryCatch(as.integer(ctrl$maxit),
                         warning = function(w) NA, error = function(e) NA)
  if (any(length(ctrl$maxit) != 1L, is.na(ctrl$maxit)))
    stop("'maxit' is not an integer")
  if (ctrl$maxit < 1)
    stop("'maxit' is not a positive integer")

  # tol: numeric (positive, finite)
  ctrl$tol <- tryCatch(as.numeric(ctrl$tol),
                          warning = function(w) NA, error = function(e) NA)
  if (any(length(ctrl$tol) != 1L, is.na(ctrl$tol),
          is.infinite(ctrl$tol)))
    stop("'tol' is not a number")
  if (ctrl$tol < .Machine$double.eps)
    stop("'tol' is not a positive number")

  # th: numeric (positive, finite)
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

check_weights <- function(weights, nw) {
  if (!is.numeric(weights))
    stop("'weights' must be a numeric vector")
  w <- as.numeric(weights)
  if (any(!is.finite(w)))
    stop("'weights' must be a finite numeric vector")
  if (any(w < 0))
    stop("negative 'weights' are not allowed")
  if (length(w) != nw)
    stop("length of 'weights' is incompatible with data")
  w <- (nw / sum(w)) * w
  w
}

check_hypothesis <- function(lhs, p) {
  if (!is.numeric(lhs) || !all(is.finite(lhs)))
    stop("'lhs' must be a finite numeric matrix")
  if (ncol(lhs) != p)
    stop("'object' and 'lhs' have incompatible dimensions")
  q <- nrow(lhs)
  if (q == 0L || q > p)
    stop("'object' and 'lhs' have incompatible dimensions")
  q
}
