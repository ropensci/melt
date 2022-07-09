#' @importFrom stats qchisq
#' @rdname confreg
#' @srrstats {G2.4, G2.4c} `as.character()` is used to `pnames`, a character
#'   vector of length two for the name of parameters.
setMethod("confreg", "EL", function(object,
                                    parm,
                                    level = 0.95,
                                    cv = NULL,
                                    npoints = 50L,
                                    control = el_control()) {
  level <- validate_level(level)
  npoints <- validate_npoints(npoints)
  stopifnot("Invalid `control` specified." = (is(control, "ControlEL")))
  npar <- getNumPar(object)
  if (npar == 0L) {
    stop("`confreg` method is not applicable to an empty model.")
  } else if (npar == 1L) {
    stop("`confreg` method is not applicable to a model with one parameter.")
  }
  est <- coef(object)
  pnames <- if (is.null(names(est))) seq(length(est)) else names(est)
  if (!missing(parm)) {
    stopifnot("Length of `parm` must be two." = (length(parm) == 2L))
    if (is.numeric(parm) && all(is.finite(parm))) {
      idx <- as.integer(parm)
      est <- est[idx]
      pnames <- pnames[idx]
      stopifnot(
        "Specify two parameters for `parm`." =
          (length(unique(est)) == 2L)
      )
    } else if (is.character(parm)) {
      stopifnot(
        "`parm` is not recognized since 'object' has no named parameters." =
          (isFALSE(is.null(names(est))))
      )
      idx <- match(parm, pnames)
      est <- est[idx]
      pnames <- pnames[idx]
      stopifnot(
        "Only one parameter specified by `parm`." =
          (length(unique(est)) == 2L)
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
      points = matrix(est, nrow = 1L), estimates = est, level = level, cv = Inf,
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
    getMethodEL(object), coef(object), getDataMatrix(object), npar, cv, idx,
    circ, control@maxit, control@maxit_l, control@tol, control@tol_l,
    control@step, control@th, control@nthreads, w
  )
  new("ConfregEL",
    points = t(circ) * cr + rep(est, each = ncol(circ)),
    estimates = est, level = level, cv = cv, pnames = as.character(pnames)
  )
})
