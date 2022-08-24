#' @rdname elmt
setMethod("elmt", "EL", function(object,
                                 rhs = NULL,
                                 lhs = NULL,
                                 alpha = 0.05,
                                 control = el_control()) {
  npar <- getNumPar(object)
  stopifnot(
    "`elmt()` method is not applicable to to an empty model." = (npar != 0L),
    "`elmt()` method is not applicable to a model with one parameter." =
      (npar != 1L),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object)))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  h <- validate_hypotheses(rhs, lhs, npar)
  qh <- head(h$q, n = length(h$q) - 1L) + 1L
  qt <- tail(h$q, n = length(h$q) - 1L)
  coefficients <- lapply(seq_along(qh), \(x) {
    drop(h$l %*% coef(object))[qh[x]:qt[x]]
  })
  alpha <- validate_alpha(alpha)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  m <- control@m
  out <- test_multiple_hypotheses(
    alpha, h$q, h$m, m, method, coef(object), getData(object), h$r, h$l, maxit,
    maxit_l, tol, tol_l, step, th, getWeights(object)
  )
  new("ELMT",
    alpha = alpha, coefficients = coefficients, statistic = out$statistic,
    df = diff(h$q), cv = out$cv, pval = out$pval, calibrate = "mvchisq"
  )
})

#' @rdname elmt
setMethod("elmt", "QGLM", function(object,
                                   rhs = NULL,
                                   lhs = NULL,
                                   alpha = 0.05,
                                   control = el_control()) {
  # `npar` includes the dispersion parameter
  p <- getNumPar(object) - 1L
  stopifnot(
    "`elmt()` method is not applicable to to an empty model." = (p != 0L),
    "`elmt()` method is not applicable to a model with one parameter." =
      (p != 1L),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object)))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  h <- validate_hypotheses(rhs, lhs, p)
  qh <- head(h$q, n = length(h$q) - 1L) + 1L
  qt <- tail(h$q, n = length(h$q) - 1L)
  coefficients <- lapply(seq_along(qh), \(x) {
    drop(h$l %*% coef(object))[qh[x]:qt[x]]
  })
  alpha <- validate_alpha(alpha)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  m <- control@m
  out <- test_multiple_hypotheses(
    alpha, h$q, h$m, m, method, c(coef(object), object@dispersion),
    getData(object), h$r, cbind(h$l, 0), maxit, maxit_l, tol, tol_l, step, th,
    getWeights(object)
  )
  new("ELMT",
    alpha = alpha, coefficients = coefficients, statistic = out$statistic,
    df = diff(h$q), cv = out$cv, pval = out$pval, calibrate = "mvchisq"
  )
})
