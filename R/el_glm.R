#' @export
el_glm <- function(formula, family = gaussian, data, weights, control = list(),
                   model = TRUE, contrasts = NULL, ...) {
  cl <- match.call()

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if (missing(data))
    data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  if (length(dim(y)) == 1L) {
    nm <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nm))
      names(y) <- nm
  }
  #
  if (is.matrix(y))
    stop("'el_glm' does not support multiple responses")
  #

  # what happens if model is empty
  if (is.empty.model(mt)) {
    # out <- list(optim = list(), npar = 0L, log.prob = numeric(),
    #             loglik = numeric(), coefficients = numeric(), df = 0L,
    #             residuals = y, fitted.values = 0 * y, na.action = action,
    #             xlevels = .getXlevels(mt, mf), call = cl, terms = mt)
    # if (keep.data)
    #   out$data.matrix <- mm
    # class(out) <- c("el_glm", "el_test")
    # return(out)
    return("empty model")
  }
  #

  x <- model.matrix(mt, mf, contrasts)


  mm <- cbind(y, x)
  intercept <- attr(mt, "intercept")
  ##
  # fitting process comes here
  aa <- glm.fit(x, y,
          # weights = weights,
          # weights = rep.int(1, nobs),
          # weights = rep.int(1, NROW(x)),
          start = NULL,
          etastart = NULL,
          mustart = NULL,
          # offset = rep.int(0, nobs),
          family = family,
          control = list(),
          intercept = attr(mt, "intercept") > 0L,
          singular.ok = TRUE)$coefficients

  ##
  optcfg <- check_control(control)
  if (missing(weights)) {
    out <- glm_("logit", mm, aa, intercept, optcfg$maxit, optcfg$tol, optcfg$th)
  } else {
    # w <- check_weights(weights, NROW(mm))
    # out <- glm_(mm, w, intercept, optcfg$maxit, optcfg$tol, optcfg$th)
    # out$weights <- w
  }

  if (model)
    out$data.matrix <- mm
  out$na.action <- attr(mf, "na.action")
  # structure(c(fit, list(call = cal, formula = formula, terms = mt,
  #                       data = data, offset = offset, control = control,
  #                       method = method,
  #                       contrasts = attr(X, "contrasts"),
  #                       xlevels = .getXlevels(mt, mf))),
  #           class = c(fit$class, c("glm", "lm")))
  out$coefficients <- setNames(out$coefficients, colnames(x))
  out$call <- cl
  out$terms <- mt
  out
  # 1
}

el_glm2 <- function(formula, family = gaussian, data,
                   weights, subset, na.action, start = NULL, etastart,
                   mustart, offset, control = list(...), model = TRUE,
                   method = "glm.fit", x = FALSE, y = TRUE, singular.ok = TRUE,
                   contrasts = NULL, ...) {
  cal <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame"))
    return(mf)
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  if (identical(method, "glm.fit"))
    control <- do.call("glm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  fit <- eval(call(if (is.function(method)) "method" else method,
                   x = X, y = Y, weights = weights, start = start,
                   etastart = etastart, mustart = mustart, offset = offset,
                   family = family, control = control,
                   intercept = attr(mt, "intercept") > 0L,
                   singular.ok = singular.ok))
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- eval(call(if (is.function(method)) "method" else method,
                      x = X[, "(Intercept)", drop = FALSE], y = Y,
                      mustart = fit$fitted.values, weights = weights,
                      offset = offset, family = family, control = control,
                      intercept = TRUE))
    if (!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  structure(c(fit, list(call = cal, formula = formula, terms = mt,
                        data = data, offset = offset, control = control,
                        method = method,
                        contrasts = attr(X, "contrasts"),
                        xlevels = .getXlevels(mt, mf))),
            class = c(fit$class, c("glm", "lm")))
}
