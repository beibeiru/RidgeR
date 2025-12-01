#' @title Secreted protein activity inference
#' @description Infer the activity of over 1000 secreted proteins from tumor gene expression profiles.
#' @param Y Gene expression matrix with gene symbol (row) x sample (column).
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomizations in the permutation test, with a default value 1000.
#' @return
#'
#' A list with four items, each is a matrix.
#' beta: regression coefficients
#' se: standard errors of coefficients
#' zscore: beta/se
#' pvalue: statistical significance
#'
#' @rdname SecAct.inference.gsl.new
#' @export
#'
SecAct.inference.gsl.new <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000)
{
  if(SigMat=="SecAct")
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }else{
    X <- read.table(SigMat,sep="\t",check.names=F)
  }

  olp <- intersect(row.names(Y),row.names(X))
  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])

  X <- scale(X)
  Y <- scale(Y)

  n <- length(olp)
  p <- ncol(X)
  m <- ncol(Y)

  res <- .C("ridge_gsl_reg",
            X=as.double(t(X)),
            Y=as.double(t(Y)),
            as.integer(n),
            as.integer(p),
            as.integer(m),
            as.double(lambda),
            as.double(nrand),
            beta=double(p*m),
            se=double(p*m),
            zscore=double(p*m),
            pvalue=double(p*m),
            df_out = double(1)
  )

  beta <- matrix(res$beta,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  se <- matrix(res$se,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  zscore <- matrix(res$zscore,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  pvalue <- matrix(res$pvalue,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))

  res <- list(beta=beta, se=se, zscore=zscore, pvalue=pvalue)

  res
}

#' @title Secreted protein activity inference
#' @description Infer the activity of over 1000 secreted proteins from tumor gene expression profiles.
#' @param Y Gene expression matrix with gene symbol (row) x sample (column).
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor in the ridge regression.
#' @param nrand Number of randomizations in the permutation test, with a default value 1000.
#' @return
#'
#' A list with four items, each is a matrix.
#' beta: regression coefficients
#' se: standard errors of coefficients
#' zscore: beta/se
#' pvalue: statistical significance
#'
#' @rdname SecAct.inference.gsl.old
#' @export
#'
SecAct.inference.gsl.old <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000)
{
  if(SigMat=="SecAct")
  {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  }else{
    X <- read.table(SigMat,sep="\t",check.names=F)
  }

  olp <- intersect(row.names(Y),row.names(X))
  X <- as.matrix(X[olp,,drop=F])
  Y <- as.matrix(Y[olp,,drop=F])

  X <- scale(X)
  Y <- scale(Y)

  n <- length(olp)
  p <- ncol(X)
  m <- ncol(Y)

  res <- .C("ridgeReg",
            X=as.double(t(X)),
            Y=as.double(t(Y)),
            as.integer(n),
            as.integer(p),
            as.integer(m),
            as.double(lambda),
            as.double(nrand),
            beta=double(p*m),
            se=double(p*m),
            zscore=double(p*m),
            pvalue=double(p*m)
  )

  beta <- matrix(res$beta,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  se <- matrix(res$se,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  zscore <- matrix(res$zscore,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))
  pvalue <- matrix(res$pvalue,byrow=T,ncol=m,dimnames=list(colnames(X),colnames(Y)))

  res <- list(beta=beta, se=se, zscore=zscore, pvalue=pvalue)

  res
}

