#' @title Secreted protein activity inference (Legacy Version)
#' @description Old single-threaded implementation.
#' @param Y Gene expression matrix.
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @export
SecAct.inference.gsl.old <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000)
{
  if(SigMat=="SecAct") {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  } else {
    X <- read.table(SigMat,sep="\t",check.names=F)
  }

  olp <- intersect(row.names(Y),row.names(X))
  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])
  X <- scale(X); Y <- scale(Y)
  X[is.na(X)] <- 0; Y[is.na(Y)] <- 0

  n <- length(olp); p <- ncol(X); m <- ncol(Y)

  res <- .C("ridgeReg",
            X=as.double(t(X)),
            Y=as.double(t(Y)),
            as.integer(n), as.integer(p), as.integer(m),
            as.double(lambda), as.double(nrand),
            beta=double(p*m), se=double(p*m), zscore=double(p*m), pvalue=double(p*m)
  )

  formatter <- function(v) matrix(v, byrow=T, ncol=m, dimnames=list(colnames(X), colnames(Y)))
  list(beta=formatter(res$beta), se=formatter(res$se), zscore=formatter(res$zscore), pvalue=formatter(res$pvalue))
}

#' @title Secreted protein activity inference (Optimized Version)
#' @description Fast multi-threaded implementation using OpenMP.
#' @param Y Gene expression matrix.
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @param ncores Number of cores (default: all available).
#' @export
SecAct.inference.gsl.new <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000, ncores=NULL)
{
  if(SigMat=="SecAct") {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  } else {
    X <- read.table(SigMat,sep="\t",check.names=F)
  }

  olp <- intersect(row.names(Y),row.names(X))
  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])
  X <- scale(X); Y <- scale(Y)
  X[is.na(X)] <- 0; Y[is.na(Y)] <- 0

  n <- length(olp); p <- ncol(X); m <- ncol(Y)
  
  if(is.null(ncores)) ncores <- parallel::detectCores(logical=FALSE)
  if(is.na(ncores)) ncores <- 1

  res <- .C("ridgeRegFast",
            X=as.double(t(X)),
            Y=as.double(t(Y)),
            as.integer(n), as.integer(p), as.integer(m),
            as.double(lambda), as.double(nrand),
            beta=double(p*m), se=double(p*m), zscore=double(p*m), pvalue=double(p*m),
            as.integer(ncores)
  )

  formatter <- function(v) matrix(v, byrow=T, ncol=m, dimnames=list(colnames(X), colnames(Y)))
  list(beta=formatter(res$beta), se=formatter(res$se), zscore=formatter(res$zscore), pvalue=formatter(res$pvalue))
}
