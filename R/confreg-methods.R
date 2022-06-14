#' @rdname confreg
setMethod(
  "confreg", "EL",
  function(object, parm, level = 0.95, cv = qchisq(level, 2L), npoints = 50L,
           control = el_control()) {
    est <- coef(object)
    if (length(est) == 0L) {
      stop("'confreg' method is not applicable to an empty model")
    } else if (length(est) == 1L) {
      stop("'confreg' method is not applicable to a model with one parameter")
    }
    pnames <- if (is.null(names(est))) seq(length(est)) else names(est)
    if (!missing(parm)) {
      if (length(parm) != 2L) {
        stop("length of 'parm' must be two")
      }
      if (is.numeric(parm) && all(is.finite(parm))) {
        idx <- as.integer(parm)
        est <- est[idx]
        pnames <- pnames[idx]
        if (length(unique(est)) != 2L) {
          stop("only one parameter specified by 'parm'")
        }
      } else if (is.character(parm)) {
        if (is.null(names(est))) {
          stop(
            "'parm' is not recognized since 'object' has no named parameters"
          )
        }
        idx <- match(parm, pnames)
        est <- est[idx]
        pnames <- pnames[idx]
        if (length(unique(est)) != 2L) {
          stop("only one parameter specified by 'parm'")
        }
      } else {
        stop("invalid 'parm' specified")
      }
    } else {
      idx <- c(1L, 2L)
      est <- est[idx]
      pnames <- pnames[idx]
    }
    if (!missing(level) &&
        (length(level) != 1L || !is.finite(level) || level < 0 || level > 1)) {
      stop("'level' must be a number between 0 and 1")
    }
    if (isTRUE(all.equal(level, 0))) {
      return(new("ConfregEL",
                 points = est, estimates = est, level = level, cv = cv,
                 pnames = c("x", "y")
      ))
    } else if (isTRUE(all.equal(level, 1))) {
      stop("'level' must be a number between 0 and 1")
    }
    if (missing(level) && !missing(cv)) {
      level <- NA_real_
    }
    cv <- check_cv(cv, control@th)
    npoints <- as.integer(npoints)
    if (npoints <= 0) {
      stop("'npoints' must be a positive integer")
    }
    w <- if (is.null(object@weights)) numeric(length = 0L) else object@weights
    ang <- seq(0, 2 * pi, length.out = npoints + 1L)[-(npoints + 1L)]
    circ <- rbind(cos(ang), sin(ang))
    if (!is(control, "ControlEL")) {
      stop("invalid 'control' specified")
    }
    cr <- confreg_(
      getMethodEL(object), coef(object), object@data, object@npar,
      cv, idx, circ, control@maxit, control@maxit_l, control@tol, control@tol_l,
      control@step, control@th, control@nthreads, w
    )
    new("ConfregEL",
        points = t(circ) * cr + rep(est, each = ncol(circ)),
        estimates = est, level = level, cv = cv, pnames = as.character(pnames)
    )
  }
)
