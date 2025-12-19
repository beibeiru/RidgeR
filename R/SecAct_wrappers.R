#' @title Secreted protein activity inference (Legacy Version)
#' @description Old single-threaded implementation using Y row permutation.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.old <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000)
{
  if(SigMat=="SecAct") {
    Xfile <- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile, sep="\t", check.names=FALSE)
  } else {
    X <- read.table(SigMat, sep="\t", check.names=FALSE)
  }
  
  olp <- intersect(row.names(Y), row.names(X))
  X <- as.matrix(X[olp, , drop=FALSE])
  Y <- as.matrix(Y[olp, , drop=FALSE])
  
  # Standardize
  X <- scale(X); Y <- scale(Y)
  X[is.na(X)] <- 0; Y[is.na(Y)] <- 0
  
  n <- length(olp)
  p <- ncol(X)
  m <- ncol(Y)
  
  res <- .C("ridgeReg",
            X = as.double(t(X)),
            Y = as.double(t(Y)),
            as.integer(n), as.integer(p), as.integer(m),
            as.double(lambda), as.double(nrand),
            beta = double(p*m), 
            se = double(p*m), 
            zscore = double(p*m), 
            pvalue = double(p*m)
  )
  
  # Note: C code outputs row-major, use byrow=TRUE
  formatter <- function(v) {
    matrix(v, byrow=TRUE, nrow=p, ncol=m, 
           dimnames=list(colnames(X), colnames(Y)))
  }
  
  list(
    beta = formatter(res$beta), 
    se = formatter(res$se), 
    zscore = formatter(res$zscore), 
    pvalue = formatter(res$pvalue)
  )
}


#' @title Secreted protein activity inference (Optimized Version - Y Permutation)
#' @description Fast multi-threaded implementation using OpenMP and .Call interface.
#'              Uses Y row permutation with blocking for cache efficiency.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @param ncores Number of cores (default: all available minus 1).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.new <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000, ncores=NULL)
{
  if(SigMat=="SecAct") {
    Xfile <- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile, sep="\t", check.names=FALSE)
  } else {
    X <- read.table(SigMat, sep="\t", check.names=FALSE)
  }
  
  olp <- intersect(row.names(Y), row.names(X))
  X <- as.matrix(X[olp, , drop=FALSE])
  Y <- as.matrix(Y[olp, , drop=FALSE])
  
  # Standardize
  X <- scale(X); X[is.na(X)] <- 0
  Y <- scale(Y); Y[is.na(Y)] <- 0
  
  if(is.null(ncores)) {
    ncores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  }
  
  p <- ncol(X)
  m <- ncol(Y)
  
  # Pass matrices directly (not transposed)
  # The C code handles the column-major to row-major mapping internally
  res <- .Call("ridgeRegFast_interface",
               X,  # (genes x proteins)
               Y,  # (genes x samples)
               as.numeric(lambda),
               as.integer(nrand),
               as.integer(ncores))
  
  # C code outputs row-major (p x m), convert to R matrix with byrow=TRUE
  formatter <- function(v) {
    matrix(v, byrow=TRUE, nrow=p, ncol=m,
           dimnames=list(colnames(X), colnames(Y)))
  }
  
  list(
    beta = formatter(res$beta), 
    se = formatter(res$se), 
    zscore = formatter(res$zscore), 
    pvalue = formatter(res$pvalue)
  )
}


#' @title Secreted protein activity inference (T Column Permutation Version)
#' @description Fast multi-threaded implementation using T column permutation.
#'              Mathematically equivalent to Y row permutation but often more
#'              efficient when number of samples (m) >> number of genes (n).
#'              
#'              The key insight: T @ Y[perm,:] == T[:,perm] @ Y
#'              
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @param ncores Number of cores (default: all available minus 1).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.new2 <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000, ncores=NULL)
{
  if(SigMat=="SecAct") {
    Xfile <- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile, sep="\t", check.names=FALSE)
  } else {
    X <- read.table(SigMat, sep="\t", check.names=FALSE)
  }
  
  olp <- intersect(row.names(Y), row.names(X))
  X <- as.matrix(X[olp, , drop=FALSE])
  Y <- as.matrix(Y[olp, , drop=FALSE])
  
  # Standardize
  X <- scale(X); X[is.na(X)] <- 0
  Y <- scale(Y); Y[is.na(Y)] <- 0
  
  if(is.null(ncores)) {
    ncores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  }
  
  p <- ncol(X)
  m <- ncol(Y)
  
  # Pass matrices directly (not transposed)
  # The C code handles the column-major to row-major mapping internally
  res <- .Call("ridgeRegTperm_interface",
               X,  # (genes x proteins)
               Y,  # (genes x samples)
               as.numeric(lambda),
               as.integer(nrand),
               as.integer(ncores))
  
  # C code outputs row-major (p x m), convert to R matrix with byrow=TRUE
  formatter <- function(v) {
    matrix(v, byrow=TRUE, nrow=p, ncol=m,
           dimnames=list(colnames(X), colnames(Y)))
  }
  
  list(
    beta = formatter(res$beta), 
    se = formatter(res$se), 
    zscore = formatter(res$zscore), 
    pvalue = formatter(res$pvalue)
  )
}


#' @title Compare all three implementations
#' @description Utility function to verify all implementations produce identical results.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @param ncores Number of cores for parallel versions.
#' @param tol Tolerance for numerical comparison (default: 1e-10).
#' @return List with results from all methods and comparison statistics.
#' @export
SecAct.compare.methods <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000, 
                                    ncores=NULL, tol=1e-10) {
  cat("Running old (single-threaded) version...\n")
  t1 <- system.time({
    res_old <- SecAct.inference.gsl.old(Y, SigMat, lambda, nrand)
  })
  cat(sprintf("  Time: %.2f seconds\n", t1["elapsed"]))
  
  cat("Running new (Y permutation, multi-threaded) version...\n")
  t2 <- system.time({
    res_new <- SecAct.inference.gsl.new(Y, SigMat, lambda, nrand, ncores)
  })
  cat(sprintf("  Time: %.2f seconds\n", t2["elapsed"]))
  
  cat("Running new2 (T permutation, multi-threaded) version...\n")
  t3 <- system.time({
    res_new2 <- SecAct.inference.gsl.new2(Y, SigMat, lambda, nrand, ncores)
  })
  cat(sprintf("  Time: %.2f seconds\n", t3["elapsed"]))
  
  # Compare results
  compare_matrices <- function(m1, m2, name) {
    max_diff <- max(abs(m1 - m2), na.rm=TRUE)
    mean_diff <- mean(abs(m1 - m2), na.rm=TRUE)
    identical_flag <- max_diff < tol
    list(max_diff=max_diff, mean_diff=mean_diff, identical=identical_flag)
  }
  
  cat("\nComparing old vs new (Y permutation):\n")
  for(stat in c("beta", "se", "zscore", "pvalue")) {
    cmp <- compare_matrices(res_old[[stat]], res_new[[stat]], stat)
    cat(sprintf("  %s: max_diff=%.2e, mean_diff=%.2e, match=%s\n", 
                stat, cmp$max_diff, cmp$mean_diff, cmp$identical))
  }
  
  cat("\nComparing old vs new2 (T permutation):\n")
  for(stat in c("beta", "se", "zscore", "pvalue")) {
    cmp <- compare_matrices(res_old[[stat]], res_new2[[stat]], stat)
    cat(sprintf("  %s: max_diff=%.2e, mean_diff=%.2e, match=%s\n", 
                stat, cmp$max_diff, cmp$mean_diff, cmp$identical))
  }
  
  cat("\nComparing new vs new2:\n")
  for(stat in c("beta", "se", "zscore", "pvalue")) {
    cmp <- compare_matrices(res_new[[stat]], res_new2[[stat]], stat)
    cat(sprintf("  %s: max_diff=%.2e, mean_diff=%.2e, match=%s\n", 
                stat, cmp$max_diff, cmp$mean_diff, cmp$identical))
  }
  
  invisible(list(
    old = res_old,
    new = res_new,
    new2 = res_new2,
    times = c(old=t1["elapsed"], new=t2["elapsed"], new2=t3["elapsed"])
  ))
}
