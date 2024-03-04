#' @rdname confreg
setMethod("confreg", "EL", function(object,
                                    parm,
                                    level = 0.95,
                                    cv = NULL,
                                    npoints = 50L,
                                    control = NULL) {
  level <- assert_number(level, lower = 0, upper = 1, finite = TRUE)
  npoints <- assert_int(npoints, lower = 1, coerce = TRUE)
  npar <- getNumPar(object)
  stopifnot(
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      isFALSE(is.null(getData(object))),
    "`confreg()` is not applicable to an empty model." = getDF(object) >= 1L,
    "`confreg()` is not applicable to a model with one parameter." =
      npar >= if (is(object, "QGLM")) 3L else 2L
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    assert_class(control, "ControlEL")
  }
  est <- coef(object)
  nm <- names(est)
  pnames <- if (is.null(nm)) seq(length(est)) else nm
  if (!missing(parm)) {
    stopifnot("Length of `parm` must be two." = length(parm) == 2L)
    if (is.numeric(parm) && all(is.finite(parm))) {
      idx <- as.integer(parm)
      est <- est[idx]
      pnames <- pnames[idx]
      stopifnot(
        "Specify two parameters for `parm`." = length(unique(est)) == 2L
      )
    } else if (is.character(parm)) {
      stopifnot(
        "`parm` is not recognized since `object` has no named parameters." =
          isFALSE(is.null(nm))
      )
      idx <- match(parm, pnames)
      stopifnot(
        "`parm` does not match the parameters in `object`." =
          isFALSE(any(is.na(idx)))
      )
      est <- est[idx]
      pnames <- pnames[idx]
      stopifnot(
        "Specify two parameters for `parm`." = length(unique(est)) == 2L
      )
    } else {
      stop("Invalid `parm` specified.")
    }
  } else {
    idx <- c(1L, 2L)
    est <- est[idx]
    pnames <- pnames[idx]
  }
  if (isTRUE(all.equal(level, 0))) {
    return(new("ConfregEL",
      .Data = matrix(est, nrow = 1L), estimates = est, level = level, cv = Inf,
      pnames = as.character(pnames)
    ))
  } else if (isTRUE(all.equal(level, 1))) {
    stop("`level` must be a number between 0 and 1.")
  }
  if (is.null(cv)) {
    cv <- qchisq(level, 2L)
  } else {
    cv <- validate_cv(cv, control@th)
    level <- NA_real_
  }
  w <- getWeights(object)
  ang <- seq(0, 2 * pi, length.out = npoints + 1L)[-(npoints + 1L)]
  circ <- rbind(cos(ang), sin(ang))
  cr <- compute_confidence_region(
    getMethodEL(object), getEstimates(object), getData(object), npar, cv, idx,
    circ, control@maxit, control@maxit_l, control@tol, control@tol_l,
    control@step, control@th, control@nthreads, w
  )
  new("ConfregEL",
    .Data = t(circ) * cr + rep(est, each = ncol(circ)),
    estimates = est, level = level, cv = cv, pnames = as.character(pnames)
  )
})
