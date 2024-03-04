#' @rdname elmt
setMethod("elmt", "EL", function(object,
                                 rhs = NULL,
                                 lhs = NULL,
                                 alpha = 0.05,
                                 control = NULL) {
  if (is(object, "QGLM")) {
    p <- getNumPar(object) - 1L
    pnames <- names(getOptim(object)$par[-getNumPar(object)])
  } else {
    p <- getNumPar(object)
    pnames <- names(getOptim(object)$par)
  }
  stopifnot(
    "`elmt()` is not applicable to to an empty model." = getDF(object) >= 1L,
    "`elmt()` is not applicable to a model with one parameter." = p != 1L,
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      isFALSE(is.null(getData(object)))
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    assert_class(control, "ControlEL")
  }
  h <- validate_hypotheses(rhs, lhs, p, pnames)
  colnames(h$l) <- pnames
  l <- if (is(object, "QGLM")) cbind(h$l, 0) else h$l
  qh <- head(h$q, n = length(h$q) - 1L) + 1L
  qt <- tail(h$q, n = length(h$q) - 1L)
  estimates <- lapply(seq_along(qh), \(x) {
    drop(h$l %*% coef(object))[qh[x]:qt[x]]
  })
  alpha <- assert_number(alpha, lower = 0, upper = 1, finite = TRUE)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  m <- control@m
  out <- test_multiple_hypotheses(
    alpha, h$q, h$m, m, method, getEstimates(object), getData(object), h$r, l,
    maxit, maxit_l, tol, tol_l, step, th, getWeights(object)
  )
  new("ELMT",
    estimates = estimates, statistic = out$statistic, df = diff(h$q),
    pval = out$pval, cv = out$cv, rhs = h$r, lhs = l, alpha = alpha,
    calibrate = "mvchisq", weights = getWeights(object),
    coefficients = getEstimates(object), method = method,
    data = if (control@keep_data) getData(object) else NULL, control = control
  )
})
