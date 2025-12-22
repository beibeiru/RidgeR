## ================================
## Internal helper functions
## ================================

.load_signature <- function(SigMat) {
  if (SigMat == "SecAct") {
    Xfile <- file.path(system.file(package = "RidgeR"), "extdata/SecAct.tsv.gz")
    if (!file.exists(Xfile)) {
      stop("Default signature matrix not found.")
    }
    read.table(Xfile, sep = "\t", check.names = FALSE)
  } else {
    read.table(SigMat, sep = "\t", check.names = FALSE)
  }
}

.prepare_matrices <- function(X, Y) {
  olp <- intersect(rownames(Y), rownames(X))
  X <- scale(as.matrix(X[olp, , drop = FALSE]))
  Y <- scale(as.matrix(Y[olp, , drop = FALSE]))

  X[is.na(X)] <- 0
  Y[is.na(Y)] <- 0

  list(X = X, Y = Y)
}

.format_result <- function(v, X, Y) {
  matrix(
    v,
    nrow = ncol(X),
    ncol = ncol(Y),
    dimnames = list(colnames(X), colnames(Y))
  )
}

.default_ncores <- function(ncores) {
  if (is.null(ncores)) {
    max(1, parallel::detectCores(logical = FALSE) - 1)
  } else {
    ncores
  }
}

## ================================
## Inference implementations
## ================================

#' @title Secreted protein activity inference (Legacy .C Version)
#' @description Original implementation using the legacy .C interface.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.legacy <- function(Y, SigMat = "SecAct",
                                       lambda = 5e+05, nrand = 1000) {
  X <- .load_signature(SigMat)
  prep <- .prepare_matrices(X, Y)

  X <- prep$X
  Y <- prep$Y

  n <- as.integer(nrow(X))
  p <- as.integer(ncol(X))
  m <- as.integer(ncol(Y))
  len <- p * m

  res <- .C(
    "ridgeReg",
    as.double(X),
    as.double(Y),
    n,
    p,
    m,
    as.double(lambda),
    as.double(nrand),
    beta   = double(len),
    se     = double(len),
    zscore = double(len),
    pvalue = double(len),
    PACKAGE = "RidgeR"
  )

  list(
    beta   = .format_result(res$beta,   X, Y),
    se     = .format_result(res$se,     X, Y),
    zscore = .format_result(res$zscore, X, Y),
    pvalue = .format_result(res$pvalue, X, Y)
  )
}

#' @title Secreted protein activity inference (Legacy Version)
#' @description Old single-threaded implementation using Y row permutation.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.old <- function(Y, SigMat = "SecAct",
                                   lambda = 5e+05, nrand = 1000) {
  X <- .load_signature(SigMat)
  prep <- .prepare_matrices(X, Y)

  res <- .Call(
    "ridgeReg_old_interface",
    prep$X,
    prep$Y,
    as.numeric(lambda),
    as.integer(nrand),
    PACKAGE = "RidgeR"
  )

  list(
    beta   = .format_result(res$beta,   prep$X, prep$Y),
    se     = .format_result(res$se,     prep$X, prep$Y),
    zscore = .format_result(res$zscore, prep$X, prep$Y),
    pvalue = .format_result(res$pvalue, prep$X, prep$Y)
  )
}

#' @title Secreted protein activity inference (Legacy Version - T Permutation)
#' @description Single-threaded implementation using cache-friendly T column permutation.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.old2 <- function(Y, SigMat = "SecAct",
                                    lambda = 5e+05, nrand = 1000) {
  X <- .load_signature(SigMat)
  prep <- .prepare_matrices(X, Y)

  res <- .Call(
    "ridgeRegTperm_old_interface",
    prep$X,
    prep$Y,
    as.numeric(lambda),
    as.integer(nrand),
    PACKAGE = "RidgeR"
  )

  list(
    beta   = .format_result(res$beta,   prep$X, prep$Y),
    se     = .format_result(res$se,     prep$X, prep$Y),
    zscore = .format_result(res$zscore, prep$X, prep$Y),
    pvalue = .format_result(res$pvalue, prep$X, prep$Y)
  )
}

#' @title Secreted protein activity inference (Optimized Version - Y Permutation)
#' @description Fast multi-threaded implementation using OpenMP and .Call interface.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @param ncores Number of cores (default: all available minus 1).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.new <- function(Y, SigMat = "SecAct",
                                   lambda = 5e+05, nrand = 1000,
                                   ncores = NULL) {
  X <- .load_signature(SigMat)
  prep <- .prepare_matrices(X, Y)
  ncores <- .default_ncores(ncores)

  res <- .Call(
    "ridgeRegFast_interface",
    prep$X,
    prep$Y,
    as.numeric(lambda),
    as.integer(nrand),
    as.integer(ncores),
    PACKAGE = "RidgeR"
  )

  list(
    beta   = .format_result(res$beta,   prep$X, prep$Y),
    se     = .format_result(res$se,     prep$X, prep$Y),
    zscore = .format_result(res$zscore, prep$X, prep$Y),
    pvalue = .format_result(res$pvalue, prep$X, prep$Y)
  )
}

#' @title Secreted protein activity inference (T Column Permutation Version)
#' @description Fast multi-threaded implementation using T column permutation.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @param ncores Number of cores (default: all available minus 1).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.new2 <- function(Y, SigMat = "SecAct",
                                    lambda = 5e+05, nrand = 1000,
                                    ncores = NULL) {
  X <- .load_signature(SigMat)
  prep <- .prepare_matrices(X, Y)
  ncores <- .default_ncores(ncores)

  res <- .Call(
    "ridgeRegTperm_interface",
    prep$X,
    prep$Y,
    as.numeric(lambda),
    as.integer(nrand),
    as.integer(ncores),
    PACKAGE = "RidgeR"
  )

  list(
    beta   = .format_result(res$beta,   prep$X, prep$Y),
    se     = .format_result(res$se,     prep$X, prep$Y),
    zscore = .format_result(res$zscore, prep$X, prep$Y),
    pvalue = .format_result(res$pvalue, prep$X, prep$Y)
  )
}

#' @title Compare all five implementations
#' @description Utility function to verify all implementations produce identical results.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @param ncores Number of cores for parallel versions.
#' @param tol Tolerance for numerical comparison (default: 1e-10).
#' @return List with results from all methods and comparison statistics.
#' @export
SecAct.compare.methods <- function(Y, SigMat = "SecAct",
                                   lambda = 5e+05, nrand = 100,
                                   ncores = NULL, tol = 1e-10) {

  time_it <- function(expr) system.time(expr)["elapsed"]

  cat("Running gsl.legacy (.C legacy)...\n")
  t0 <- time_it(res_leg  <- SecAct.inference.gsl.legacy(Y, SigMat, lambda, nrand))

  cat("Running gsl.old (single-threaded, Y row permutation)...\n")
  t1 <- time_it(res_old  <- SecAct.inference.gsl.old(Y, SigMat, lambda, nrand))

  cat("Running gsl.old2 (single-threaded, T column permutation)...\n")
  t2 <- time_it(res_old2 <- SecAct.inference.gsl.old2(Y, SigMat, lambda, nrand))

  cat("Running gsl.new (multi-threaded, Y row permutation)...\n")
  t3 <- time_it(res_new  <- SecAct.inference.gsl.new(Y, SigMat, lambda, nrand, ncores))

  cat("Running gsl.new2 (multi-threaded, T column permutation)...\n")
  t4 <- time_it(res_new2 <- SecAct.inference.gsl.new2(Y, SigMat, lambda, nrand, ncores))

  equal <- function(a, b) max(abs(a - b), na.rm = TRUE) < tol

  cat("\n=== Consistency Check (z-score) ===\n")
  cat(sprintf("legacy vs old  : %s\n", equal(res_leg$zscore,  res_old$zscore)))
  cat(sprintf("old    vs new  : %s\n", equal(res_old$zscore,  res_new$zscore)))
  cat(sprintf("old2   vs new2 : %s\n", equal(res_old2$zscore, res_new2$zscore)))
  cat(sprintf("new    vs new2 : %s\n", equal(res_new$zscore,  res_new2$zscore)))

  cat("\n=== Timing Summary ===\n")
  cat(sprintf("  gsl.legacy : %.2fs\n", t0))
  cat(sprintf("  gsl.old    : %.2fs\n", t1))
  cat(sprintf("  gsl.old2   : %.2fs\n", t2))
  cat(sprintf("  gsl.new    : %.2fs\n", t3))
  cat(sprintf("  gsl.new2   : %.2fs\n", t4))

  invisible(list(
    legacy = res_leg,
    old    = res_old,
    old2   = res_old2,
    new    = res_new,
    new2   = res_new2
  ))
}
