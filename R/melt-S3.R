#' @export
print.el_aov <- function(x, ...)
{
  stopifnot(inherits(x, "melt"))
  cat("Call:\n")
  dput(x$call, control = NULL)
  cat("\nminimizer:\n")
  cat(format(round(x$optim$par, 4), scientific = FALSE))
  cat("\n\n")
  cat("statistic:\n")
  cat(format(round(x$optim$n2logLR, 4), scientific = FALSE))
  cat("\n\n")
}

# print.el_test <- function(x, ...) {
#   stopifnot(inherits(x, "melt"))
#   cat("\n")
#   cat("Empirical Likelihood Hypothesis Testing\n\n")
#   cat("minimizer:\n")
#   cat(format(round(x$optim$par, 4), scientific = F))
#   cat("\n\n")
#   cat("statistic:\n")
#   cat(format(round(x$optim$n2logLR, 4), scientific = F))
#   cat("\n\n")
# }

#' @importFrom stats printCoefmat
#' @export
print.pairwise <- function(x, ...) {
  stopifnot(inherits(x, "melt"))
  cat("\n")
  cat("Empirical Likelihood Multiple Hypothesis Testing\n\n")
  # set row names
  if (is.null(x$control)) {
    cat("Test: all pairwise comparisons\n\n")
    rname <- vector("character", length = 0)
    for (i in 1L:(length(x$trt) - 1L)) {
      for (j in (i + 1L):length(x$trt)) {
        rname <- c(rname, paste(x$trt[i], "-", x$trt[j]))
      }
    }
  } else {
    cat("Test: comparisons with control\n\n")
    diff <- setdiff(x$trt, x$control)
    rname <- vector("character", length = length(diff))
    for (i in 1L:length(diff)) {
      rname[i] <- paste(diff[i], "-", x$control)
    }
  }
  out <- data.frame(row.names = rname, estimate  = x$estimate,
                    statistic = x$statistic, lwr.ci = x$lower,
                    upr.ci = x$upper,
                    p.adj = round(x$p.adj, 4))
  printCoefmat(out, digits = min(4L, getOption("digits")), cs.ind = c(1, 3, 4),
               tst.ind = 2L, dig.tst = min(3L, getOption("digits")),
               P.values = TRUE, has.Pvalue = TRUE, eps.Pvalue = 1e-03)
  cat("\n")
  cat(paste(c("k", "level", "method", "cutoff"),
            c(x$k, x$level, x$method, round(x$cutoff, 4)),
            collapse = ", ", sep = ": "))
  cat("\n\n")
}

# check_hypotheses <- function(hypotheses, rhs, trt) {
#   # hypotheses must be either character vector or list of numeric matrices
#   match.arg(class(hypotheses), c("character", "list"))
#   # test for pairwise comparisons
#   test <- NULL
#   if (is.character(hypotheses)) {
#     if (length(hypotheses) == 1 && hypotheses == "pairwise") {
#       test <- 0
#       # all pairwise
#     } else if (length(hypotheses) == 2 && hypotheses[1] == "pairwise") {
#       if (hypotheses[2] %in% trt) {
#         test <- which(trt == hypotheses[2])
#       } else {
#         stop(paste("there is no treatment named", hypotheses[2]))
#       }
#     } else {
#       stop("invalid hypotheses specified")
#     }
#     return(test)
#   }
#
#   # check validity of hypotheses and rhs
#   if (!all(sapply(hypotheses, function(x) is.numeric(x) && is.matrix(x)))) {
#     stop("hypotheses must consists of numeric matrices")
#   } else if (!all(sapply(hypotheses, function(x) {
#     nrow(x) <= ncol(x) && ncol(x) == length(trt)}))) {
#     stop("invalid matrix dimensions for hypotheses")
#   }
#
#   if (!is.null(rhs)) {
#     if (!is.list(rhs)) {
#       stop("invalid rhs specified")
#     } else if (length(hypotheses) != length(rhs)) {
#       stop("number of hypotheses and rhs do not match")
#     } else if (!all(sapply(rhs, function(x) is.numeric(x) && is.vector(x)))) {
#       stop("rhs must consists of numeric vectors")
#     }
#   }
#   test
# }

# elmulttest_rownames <- function(m, trt, test) {
#   # if not pairwise comparison, return general hypotheses names
#   if (is.null(test)) {
#     return(paste0("H_", 1:m))
#   }
#
#   # pairwise comparisons
#   if (length(test) == 2) {
#     # comparisons with control
#     diff <- setdiff(trt, test[2])
#     row_names <- vector("character", length = length(diff))
#     for (i in 1:length(diff)) {
#       row_names[i] <- paste(diff[i], "-", test[2])
#     }
#     return(row_names)
#   } else {
#     # all pairwise comparisons
#     row_names <- vector("character", length = 0)
#     for (i in 1:(length(trt) - 1)) {
#       for (j in (i + 1):length(trt)) {
#         row_names <- c(row_names, paste(trt[i], "-", trt[j]))
#       }
#     }
#     row_names
#   }
# }

