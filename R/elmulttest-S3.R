check_hypotheses <- function(hypotheses, rhs, trt) {
  # hypotheses must be either character vector or list of numeric matrices
  match.arg(class(hypotheses), c("character", "list"))
  # test for pairwise comparisons
  test <- NULL
  if (is.character(hypotheses)) {
    if (length(hypotheses) == 1 && hypotheses == "pairwise") {
      test <- 0
      # all pairwise
    } else if (length(hypotheses) == 2 && hypotheses[1] == "pairwise") {
      if (hypotheses[2] %in% trt) {
        test <- which(trt == hypotheses[2])
      } else {
        stop(paste("there is no treatment named", hypotheses[2]))
      }
    } else {
      stop("invalid hypotheses specified")
    }
    return(test)
  }

  # check validity of hypotheses and rhs
  if (!all(sapply(hypotheses, function(x) is.numeric(x) && is.matrix(x)))) {
    stop("hypotheses must consists of numeric matrices")
  } else if (!all(sapply(hypotheses, function(x) {
    nrow(x) <= ncol(x) && ncol(x) == length(trt)}))) {
    stop("invalid matrix dimensions for hypotheses")
  }

  if (!is.null(rhs)) {
    if (!is.list(rhs)) {
      stop("invalid rhs specified")
    } else if (length(hypotheses) != length(rhs)) {
      stop("number of hypotheses and rhs do not match")
    } else if (!all(sapply(rhs, function(x) is.numeric(x) && is.vector(x)))) {
      stop("rhs must consists of numeric vectors")
    }
  }
  test
}

elmulttest_rownames <- function(m, trt, test) {
  # if not pairwise comparison, return general hypotheses names
  if (is.null(test)) {
    return(paste0("H_", 1:m))
  }

  # pairwise comparisons
  if (length(test) == 2) {
    # comparisons with control
    diff <- setdiff(trt, test[2])
    row_names <- vector("character", length = length(diff))
    for (i in 1:length(diff)) {
      row_names[i] <- paste(diff[i], "-", test[2])
    }
    return(row_names)
  } else {
    # all pairwise comparisons
    row_names <- vector("character", length = 0)
    for (i in 1:(length(trt) - 1)) {
      for (j in (i + 1):length(trt)) {
        row_names <- c(row_names, paste(trt[i], "-", trt[j]))
      }
    }
    row_names
  }
}

#' @export
print.pairwise <- function(x, ...) {
  stopifnot(inherits(x, "elmulttest"))
  cat("Empirical Likelihood Multiple Hypothesis Testing\n\n")
  if (all(x$test == "pairwise")) {
    cat("Test: all pairwise comparisons\n\n")
  } else {
    cat("Test: comparisons with control\n\n")
  }
  out <- data.frame(row.names =
                      elmulttest_rownames(length(x$statistic), x$trt, x$test))
  out$estimate  <- x$estimate
  out$statistic <- x$statistic
  out$lwr.ci    <- x$lower
  out$upr.ci    <- x$upper
  out$p.adj     <- x$p.adj
  print(format(round(out, 4), digits = 4))
  cat("---\n", "k: ", x$k, ", level: ", x$level, ", method: ", x$method,
      ", cutoff: ", round(x$cutoff, 4), "\n", sep = "")
}
