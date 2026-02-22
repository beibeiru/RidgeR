#' @keywords internal
group_signatures <- function(X, cor_threshold = 0.9) {
  dis <- as.dist(1 - cor(X, method = "pearson"))
  hc <- hclust(dis, method = "complete")
  group_labels <- cutree(hc, h = 1 - cor_threshold)
  newsig <- data.frame(row.names = rownames(X))
  for (j in unique(group_labels)) {
    geneGroups <- names(group_labels)[group_labels == j]
    newsig[, paste0(geneGroups, collapse = "|")] <- rowMeans(X[, geneGroups, drop = FALSE])
  }
  as.data.frame(newsig)
}

#' @keywords internal
expand_rows <- function(mat) {
  new_rows <- lapply(1:nrow(mat), function(i) {
    names <- strsplit(rownames(mat)[i], "\\|")[[1]]
    `rownames<-`(do.call(rbind, replicate(length(names), mat[i, , drop = FALSE], simplify = FALSE)), names)
  })
  do.call(rbind, new_rows)
}

#' @title Secreted protein activity inference (Legacy Version)
#' @description Old single-threaded implementation.
#' @param Y Gene expression matrix.
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @param rng_method RNG method: "srand" (default, platform-dependent C stdlib rand)
#'   or "gsl" (cross-platform GSL MT19937, matches SecActPy GSLRNG).
#' @param is.group.sig Logical; group correlated signatures before regression (default TRUE).
#' @param is.group.cor Correlation threshold for grouping (default 0.9).
#' @export
SecAct.inference.gsl.old <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000, rng_method="srand",
                                     is.group.sig=TRUE, is.group.cor=0.9)
{
  if(SigMat=="SecAct") {
    Xfile<- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile,sep="\t",check.names=F)
  } else {
    X <- read.table(SigMat,sep="\t",check.names=F)
  }

  if (is.group.sig) {
    X <- group_signatures(X, cor_threshold = is.group.cor)
  }

  rng_int <- match.arg(rng_method, c("srand", "gsl"))
  rng_int <- ifelse(rng_int == "gsl", 1L, 0L)

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
            beta=double(p*m), se=double(p*m), zscore=double(p*m), pvalue=double(p*m),
            rng_method=as.integer(rng_int)
  )

  formatter <- function(v) matrix(v, byrow=T, ncol=m, dimnames=list(colnames(X), colnames(Y)))
  result <- list(beta=formatter(res$beta), se=formatter(res$se), zscore=formatter(res$zscore), pvalue=formatter(res$pvalue))

  if (is.group.sig) {
    result$beta   <- expand_rows(result$beta)
    result$se     <- expand_rows(result$se)
    result$zscore <- expand_rows(result$zscore)
    result$pvalue <- expand_rows(result$pvalue)
    ord <- sort(rownames(result$beta))
    result$beta   <- result$beta[ord, , drop = FALSE]
    result$se     <- result$se[ord, , drop = FALSE]
    result$zscore <- result$zscore[ord, , drop = FALSE]
    result$pvalue <- result$pvalue[ord, , drop = FALSE]
  }

  result
}

#' @title Secreted protein activity inference (Optimized Version)
#' @description Dispatches to one of 4 optimized variants based on the \code{method} parameter.
#'   Default is \code{"Tcol.mt"} (T-column permutation, multi-threaded).
#' @param Y Gene expression matrix.
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @param ncores Number of cores (default: all available; multi-threaded variants only).
#' @param rng_method RNG method: "srand" (default, platform-dependent C stdlib rand)
#'   or "gsl" (cross-platform GSL MT19937, matches SecActPy GSLRNG).
#' @param method Variant to use: "Tcol.mt" (default), "Tcol.st", "Yrow.mt", or "Yrow.st".
#' @param is.group.sig Logical; group correlated signatures before regression (default TRUE).
#' @param is.group.cor Correlation threshold for grouping (default 0.9).
#' @export
SecAct.inference.gsl.new <- function(Y, SigMat="SecAct", lambda=5e+05, nrand=1000, ncores=NULL, rng_method="srand",
                                     method="Tcol.mt", is.group.sig=TRUE, is.group.cor=0.9)
{
  method <- match.arg(method, c("Tcol.mt", "Tcol.st", "Yrow.mt", "Yrow.st"))

  switch(method,
    Tcol.mt = SecAct.inference.Tcol.mt(Y, SigMat=SigMat, lambda=lambda, nrand=nrand,
                                       ncores=ncores, rng_method=rng_method,
                                       is.group.sig=is.group.sig, is.group.cor=is.group.cor),
    Tcol.st = SecAct.inference.Tcol.st(Y, SigMat=SigMat, lambda=lambda, nrand=nrand,
                                       rng_method=rng_method,
                                       is.group.sig=is.group.sig, is.group.cor=is.group.cor),
    Yrow.mt = SecAct.inference.Yrow.mt(Y, SigMat=SigMat, lambda=lambda, nrand=nrand,
                                       ncores=ncores, rng_method=rng_method,
                                       is.group.sig=is.group.sig, is.group.cor=is.group.cor),
    Yrow.st = SecAct.inference.Yrow.st(Y, SigMat=SigMat, lambda=lambda, nrand=nrand,
                                       rng_method=rng_method,
                                       is.group.sig=is.group.sig, is.group.cor=is.group.cor)
  )
}


# =========================================================
# Shared helpers for the 5 new variants
# =========================================================

#' @keywords internal
.secact_preprocess <- function(Y, SigMat, is.group.sig, is.group.cor, rng_method) {
  if (SigMat == "SecAct") {
    Xfile <- file.path(system.file(package = "SecAct"), "extdata/SecAct.tsv.gz")
    X <- read.table(Xfile, sep = "\t", check.names = FALSE)
  } else {
    X <- read.table(SigMat, sep = "\t", check.names = FALSE)
  }

  if (is.group.sig) {
    X <- group_signatures(X, cor_threshold = is.group.cor)
  }

  rng_int <- match.arg(rng_method, c("srand", "gsl"))
  rng_int <- ifelse(rng_int == "gsl", 1L, 0L)

  olp <- intersect(row.names(Y), row.names(X))
  X <- as.matrix(X[olp, , drop = FALSE])
  Y <- as.matrix(Y[olp, , drop = FALSE])

  X <- scale(X); X[is.na(X)] <- 0
  Y <- scale(Y); Y[is.na(Y)] <- 0

  list(X = X, Y = Y, rng_int = rng_int)
}

#' @keywords internal
.secact_postprocess <- function(result, is.group.sig) {
  if (is.group.sig) {
    result$beta   <- expand_rows(result$beta)
    result$se     <- expand_rows(result$se)
    result$zscore <- expand_rows(result$zscore)
    result$pvalue <- expand_rows(result$pvalue)
    ord <- sort(rownames(result$beta))
    result$beta   <- result$beta[ord, , drop = FALSE]
    result$se     <- result$se[ord, , drop = FALSE]
    result$zscore <- result$zscore[ord, , drop = FALSE]
    result$pvalue <- result$pvalue[ord, , drop = FALSE]
  }
  result
}

# FIX: Use byrow=TRUE to match the C backend's row-major output layout.
# The C core writes flat output as [r * m + s] (row-major, protein x sample),
# so we must read it back row-by-row. The old dim(v) <- c(p, m) approach
# was column-major and caused proteins and samples to be scrambled.
#' @keywords internal
.secact_format_call <- function(res, X, Y) {
  p <- ncol(X)
  m <- ncol(Y)
  formatter <- function(v) {
    matrix(v, nrow = p, ncol = m, byrow = TRUE,
           dimnames = list(colnames(X), colnames(Y)))
  }
  list(
    beta   = formatter(res$beta),
    se     = formatter(res$se),
    zscore = formatter(res$zscore),
    pvalue = formatter(res$pvalue)
  )
}


# =========================================================
# 5 New Variants
# =========================================================

#' @title Y-row permutation, single-threaded (.Call)
#' @description Single-threaded Y-row permutation using .Call interface (64-bit indexing).
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param rng_method RNG method: "srand" or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.Yrow.st <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     rng_method = "srand", is.group.sig = TRUE, is.group.cor = 0.9)
{
  pre <- .secact_preprocess(Y, SigMat, is.group.sig, is.group.cor, rng_method)

  res <- .Call("ridgeRegFast_interface",
               pre$X, pre$Y,
               as.numeric(lambda),
               as.integer(nrand),
               1L,  # ncores = 1 (single-threaded)
               as.integer(pre$rng_int))

  result <- .secact_format_call(res, pre$X, pre$Y)
  .secact_postprocess(result, is.group.sig)
}

#' @title Y-row permutation, multi-threaded (.Call)
#' @description Multi-threaded Y-row permutation using .Call interface with OpenMP.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param ncores Number of CPU cores (NULL = auto-detect).
#' @param rng_method RNG method: "srand" or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.Yrow.mt <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     ncores = NULL, rng_method = "srand",
                                     is.group.sig = TRUE, is.group.cor = 0.9)
{
  pre <- .secact_preprocess(Y, SigMat, is.group.sig, is.group.cor, rng_method)

  if (is.null(ncores)) {
    ncores <- parallel::detectCores(logical = FALSE)
  }

  res <- .Call("ridgeRegFast_interface",
               pre$X, pre$Y,
               as.numeric(lambda),
               as.integer(nrand),
               as.integer(ncores),
               as.integer(pre$rng_int))

  result <- .secact_format_call(res, pre$X, pre$Y)
  .secact_postprocess(result, is.group.sig)
}

#' @title T-column permutation, single-threaded (.Call)
#' @description Single-threaded T-column permutation using .Call interface.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param rng_method RNG method: "srand" or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.Tcol.st <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     rng_method = "srand", is.group.sig = TRUE, is.group.cor = 0.9)
{
  pre <- .secact_preprocess(Y, SigMat, is.group.sig, is.group.cor, rng_method)

  res <- .Call("ridgeRegFastTcol_interface",
               pre$X, pre$Y,
               as.numeric(lambda),
               as.integer(nrand),
               1L,  # ncores = 1 (single-threaded)
               as.integer(pre$rng_int))

  result <- .secact_format_call(res, pre$X, pre$Y)
  .secact_postprocess(result, is.group.sig)
}

#' @title T-column permutation, multi-threaded (.Call)
#' @description Multi-threaded T-column permutation using .Call interface with OpenMP.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param ncores Number of CPU cores (NULL = auto-detect).
#' @param rng_method RNG method: "srand" or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.Tcol.mt <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     ncores = NULL, rng_method = "srand",
                                     is.group.sig = TRUE, is.group.cor = 0.9)
{
  pre <- .secact_preprocess(Y, SigMat, is.group.sig, is.group.cor, rng_method)

  if (is.null(ncores)) {
    ncores <- parallel::detectCores(logical = FALSE)
  }

  res <- .Call("ridgeRegFastTcol_interface",
               pre$X, pre$Y,
               as.numeric(lambda),
               as.integer(nrand),
               as.integer(ncores),
               as.integer(pre$rng_int))

  result <- .secact_format_call(res, pre$X, pre$Y)
  .secact_postprocess(result, is.group.sig)
}

#' @title Pure R naive implementation with C-generated permutation table
#' @description Pure R linear algebra with C-generated permutation table for
#'   identical RNG sequences. Uses T-column permutation approach.
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param rng_method RNG method: "srand" or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.naive <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                   rng_method = "srand", is.group.sig = TRUE, is.group.cor = 0.9)
{
  pre <- .secact_preprocess(Y, SigMat, is.group.sig, is.group.cor, rng_method)
  X <- pre$X
  Y <- pre$Y

  n <- nrow(X)
  p <- ncol(X)
  m <- ncol(Y)

  # 1. Compute projection matrix: T = (X'X + lambda*I)^-1 X'
  T_mat <- solve(crossprod(X) + lambda * diag(p), t(X))  # (p x n)

  # 2. Observed beta
  beta <- T_mat %*% Y  # (p x m)

  # 3. Get permutation table from C (same RNG sequence)
  perm_table <- .Call("generate_perm_table",
                      as.integer(n),
                      as.integer(nrand),
                      as.integer(pre$rng_int))
  # perm_table is (nrand x n) integer matrix, 0-indexed

  # 4. Permutation loop (T-column approach in pure R)
  aver <- matrix(0.0, nrow = p, ncol = m)
  aver_sq <- matrix(0.0, nrow = p, ncol = m)
  pvalue_counts <- matrix(0.0, nrow = p, ncol = m)
  abs_beta <- abs(beta)

  for (i in seq_len(nrand)) {
    # Forward permutation from table (0-indexed -> 1-indexed)
    fwd_perm <- perm_table[i, ] + 1L

    # Compute inverse permutation: inv_perm[fwd_perm[j]] = j
    inv_perm <- integer(n)
    inv_perm[fwd_perm] <- seq_len(n)

    # T-column permutation: permute columns of T
    beta_rand <- T_mat[, inv_perm, drop = FALSE] %*% Y

    pvalue_counts <- pvalue_counts + (abs(beta_rand) >= abs_beta) * 1.0
    aver <- aver + beta_rand
    aver_sq <- aver_sq + beta_rand^2
  }

  # 5. Finalize statistics (same formula as C)
  mean_rand <- aver / nrand
  var_rand <- (aver_sq / nrand) - mean_rand^2
  var_rand[var_rand < 0] <- 0
  se <- sqrt(var_rand)
  zscore <- ifelse(se > 1e-12, (beta - mean_rand) / se, 0.0)
  pvalue <- (pvalue_counts + 1.0) / (nrand + 1.0)

  # Format output
  rownames(beta) <- colnames(X); colnames(beta) <- colnames(Y)
  rownames(se) <- colnames(X); colnames(se) <- colnames(Y)
  rownames(zscore) <- colnames(X); colnames(zscore) <- colnames(Y)
  rownames(pvalue) <- colnames(X); colnames(pvalue) <- colnames(Y)

  result <- list(beta = beta, se = se, zscore = zscore, pvalue = pvalue)
  .secact_postprocess(result, is.group.sig)
}
