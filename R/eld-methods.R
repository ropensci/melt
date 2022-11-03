#' @rdname eld
setMethod("eld", "EL", function(object, control = NULL) {
  stopifnot(
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object))))
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    stopifnot("Invalid `control` specified." = is(control, "ControlEL"))
  }
  new("ELD", .Data = compute_ELD(
    getMethodEL(object), coef(object), getData(object), control@maxit_l,
    control@tol_l, control@th, control@nthreads, getWeights(object)
  ))
})

#' @rdname eld
setMethod("eld", "GLM", function(object, control = NULL) {
  stopifnot(
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object))))
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    stopifnot("Invalid `control` specified." = is(control, "ControlEL"))
  }
  mm <- getData(object)
  n <- nobs(object)
  x <- mm[, -c(1L, 2L)]
  if (is.null(dim(x))) {
    dim(x) <- c(n, ncol(mm) - 2L)
  }
  y <- mm[, 2L]
  w <- weights(object)
  glm_control <- object@misc$control
  intercept <- object@misc$intercept
  offset <- object@misc$offset
  new("ELD", .Data = vapply(seq_len(n), function(i) {
    fit <- suppressWarnings(glm.fit(
      x = x[-i, ], y = y[-i], weights = w[-i], offset = offset,
      family = object@family, control = glm_control, intercept = intercept,
      singular.ok = FALSE
    ))
    -2 * n * logLR(elt(object, rhs = fit$coefficients, control = control))
  }, numeric(1L)))
})

#' @rdname eld
#' @usage NULL
setMethod("eld", "missing", function(object, control = NULL) {
  if (is.null(control)) {
    control <- el_control()
  } else {
    stopifnot("Invalid `control` specified." = is(control, "ControlEL"))
  }
  NULL
})
