#' #' Fits a generalized linear model with empirical likelihood
#' #'
#' #' Fits a generalized linear model with empirical likelihood.
#' #'
#' #' @param formula A formula object.
#' #' @param family dd
#' #' @param data A data frame containing the variables in the formula.
#' #' @param weights An optional numeric vector of weights to be used in the
#' #'   fitting process. If not provided, identical weights are applied. Otherwise,
#' #'   weighted empirical likelihood is computed.
#' #' @param control A list of control parameters. See ‘Details’ in
#' #'   \code{\link{el_eval}}.
#' #' @param model A logical. If \code{TRUE} the model matrix used for fitting is
#' #'   returned.
#' #' @param contrasts dd
#' #' @param ... For glm: arguments to be used to form the default control argument
#' #'   if it is not supplied directly.
#' #' @return A list with class \code{c("el_lm", "el_test")}.
#' #' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.”
#' #'   The Annals of Statistics 19 (4).
#' #'   \doi{10.1214/aos/1176348368}.
#' #' @seealso \link{el_eval}, \link{lht}
#' #' @examples
#' #' fit <- el_lm(mpg ~ wt, mtcars)
#' #' summary(fit)
#' #' @importFrom stats gaussian glm.fit
#' el_glm <- function(formula, family = gaussian, data, weights, control = list(),
#'                    model = TRUE, contrasts = NULL, ...) {
#'   cl <- match.call()
#'
#'   if (is.character(family))
#'     family <- get(family, mode = "function", envir = parent.frame())
#'   if (is.function(family))
#'     family <- family()
#'   if (is.null(family$family)) {
#'     print(family)
#'     stop("'family' not recognized")
#'   }
#'
#'   if (missing(data))
#'     data <- environment(formula)
#'
#'   mf <- match.call(expand.dots = FALSE)
#'   m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
#'   mf <- mf[c(1L, m)]
#'   mf$drop.unused.levels <- TRUE
#'   mf[[1L]] <- quote(stats::model.frame)
#'   mf <- eval(mf, parent.frame())
#'   mt <- attr(mf, "terms")
#'   y <- model.response(mf, "any")
#'   if (length(dim(y)) == 1L) {
#'     nm <- rownames(y)
#'     dim(y) <- NULL
#'     if (!is.null(nm))
#'       names(y) <- nm
#'   }
#'   #
#'   if (is.matrix(y))
#'     stop("'el_glm' does not support multiple responses")
#'   #
#'
#'   # what happens if model is empty
#'   if (is.empty.model(mt)) {
#'     # out <- list(optim = list(), npar = 0L, log.prob = numeric(),
#'     #             loglik = numeric(), coefficients = numeric(), df = 0L,
#'     #             residuals = y, fitted.values = 0 * y, na.action = action,
#'     #             xlevels = .getXlevels(mt, mf), call = cl, terms = mt)
#'     # if (keep.data)
#'     #   out$data.matrix <- mm
#'     # class(out) <- c("el_glm", "el_test")
#'     # return(out)
#'     return("empty model")
#'   }
#'   #
#'
#'   x <- model.matrix(mt, mf, contrasts)
#'
#'
#'   mm <- cbind(y, x)
#'   intercept <- attr(mt, "intercept")
#'   ##
#'   # fitting process comes here
#'   aa <- glm.fit(x, y,
#'           # weights = weights,
#'           # weights = rep.int(1, nobs),
#'           # weights = rep.int(1, NROW(x)),
#'           start = NULL,
#'           etastart = NULL,
#'           mustart = NULL,
#'           # offset = rep.int(0, nobs),
#'           family = family,
#'           control = list(),
#'           intercept = attr(mt, "intercept") > 0L,
#'           singular.ok = TRUE)$coefficients
#'
#'   ##
#'   optcfg <- check_control(control)
#'   if (missing(weights)) {
#'     out <- glm_("logit", mm, aa, intercept, optcfg$maxit, optcfg$tol, optcfg$th)
#'   } else {
#'     # w <- check_weights(weights, NROW(mm))
#'     # out <- glm_w_("logit", mm, w, intercept, optcfg$maxit, optcfg$tol, optcfg$th)
#'     # out$weights <- w
#'   }
#'
#'   if (model)
#'     out$data.matrix <- mm
#'   out$na.action <- attr(mf, "na.action")
#'   # structure(c(fit, list(call = cal, formula = formula, terms = mt,
#'   #                       data = data, offset = offset, control = control,
#'   #                       method = method,
#'   #                       contrasts = attr(X, "contrasts"),
#'   #                       xlevels = .getXlevels(mt, mf))),
#'   #           class = c(fit$class, c("glm", "lm")))
#'   out$coefficients <- setNames(out$coefficients, colnames(x))
#'   out$call <- cl
#'   out$terms <- mt
#'   out
#' }
#'
#'
