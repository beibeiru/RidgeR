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
  
  p <- ncol(X)
  m <- ncol(Y)
  
  res <- .Call("ridgeReg_old_interface", 
               X, Y, 
               as.numeric(lambda), 
               as.integer(nrand))
  
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


#' @title Secreted protein activity inference (Legacy Version - T Permutation)
#' @description Single-threaded implementation using cache-friendly T column permutation.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor (default: 5e+05).
#' @param nrand Number of randomizations (default: 1000).
#' @return List with beta, se, zscore, pvalue matrices (proteins x samples).
#' @export
SecAct.inference.gsl.old2 <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000)
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
  
  p <- ncol(X)
  m <- ncol(Y)
  
  res <- .Call("ridgeRegTperm_old_interface", 
               X, Y, 
               as.numeric(lambda), 
               as.integer(nrand))
  
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
  
  res <- .Call("ridgeRegFast_interface",
               X, Y,
               as.numeric(lambda),
               as.integer(nrand),
               as.integer(ncores))
  
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
  
  res <- .Call("ridgeRegTperm_interface",
               X, Y,
               as.numeric(lambda),
               as.integer(nrand),
               as.integer(ncores))
  
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

#' @title Compare all four implementations
#' @description Utility function to verify all implementations produce identical results.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Secreted protein signature matrix or "SecAct" for default.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @param ncores Number of cores for parallel versions.
#' @param tol Tolerance for numerical comparison (default: 1e-10).
#' @return List with results from all methods and comparison statistics.
#' @export
SecAct.compare.methods <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=100, 
                                    ncores=NULL, tol=1e-10) {
  cat("Running gsl.old (single-threaded, Y row permutation)...\n")
  t1 <- system.time({ res_old <- SecAct.inference.gsl.old(Y, SigMat, lambda, nrand) })
  
  cat("Running gsl.old2 (single-threaded, T column permutation)...\n")
  t2 <- system.time({ res_old2 <- SecAct.inference.gsl.old2(Y, SigMat, lambda, nrand) })
  
  cat("Running gsl.new (multi-threaded, Y row permutation)...\n")
  t3 <- system.time({ res_new <- SecAct.inference.gsl.new(Y, SigMat, lambda, nrand, ncores) })
  
  cat("Running gsl.new2 (multi-threaded, T column permutation)...\n")
  t4 <- system.time({ res_new2 <- SecAct.inference.gsl.new2(Y, SigMat, lambda, nrand, ncores) })
  
  compare_matrices <- function(m1, m2) max(abs(m1 - m2), na.rm=TRUE) < tol
  
  cat("\n=== Consistency Check ===\n")
  cat(sprintf("old  vs new : %s\n", compare_matrices(res_old$zscore, res_new$zscore)))
  cat(sprintf("old2 vs new2: %s\n", compare_matrices(res_old2$zscore, res_new2$zscore)))
  cat(sprintf("new  vs new2: %s\n", compare_matrices(res_new$zscore, res_new2$zscore)))
  
  cat("\n=== Timing Summary ===\n")
  cat(sprintf("  gsl.old  : %.2fs\n", t1["elapsed"]))
  cat(sprintf("  gsl.old2 : %.2fs\n", t2["elapsed"]))
  cat(sprintf("  gsl.new  : %.2fs\n", t3["elapsed"]))
  cat(sprintf("  gsl.new2 : %.2fs\n", t4["elapsed"]))
  
  invisible(list(old=res_old, old2=res_old2, new=res_new, new2=res_new2))
}
