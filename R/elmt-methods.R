#' @importFrom methods is
#' @importFrom stats pchisq
#' @rdname elmt
setMethod("elmt", "EL", function(object,
                                 rhs = NULL,
                                 lhs = NULL,
                                 alpha = 0.05,
                                 control = el_control()) {
  stopifnot(
    "'object' has no 'data'; fit the model with 'model' = TRUE" =
      (length(getDataMatrix(object)) > 1L),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  p <- object@npar
  h <- validate_hypotheses(rhs, lhs, p)
  # return(h)

  alpha <- validate_alpha(alpha)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th


  out <- elmt_statistic_(
    h$q, h$m, method, coef(object), getDataMatrix(object), h$r, h$l,
    maxit, maxit_l, tol, tol_l, step, th, getWeights(object)
  )
  new("MELT",
      statistic = out
  )
})
