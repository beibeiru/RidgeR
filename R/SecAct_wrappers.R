#' @keywords internal
group_signatures <- function(X, cor_threshold = 0.9) {
  dis <- as.dist(1 - cor(X, method = "pearson"))
  hc <- hclust(dis, method = "complete")
  group_labels <- cutree(hc, h = 1 - cor_threshold)
  split_names <- split(names(group_labels), group_labels)
  col_names <- vapply(split_names, paste0, character(1), collapse = "|")
  newsig <- matrix(NA_real_, nrow = nrow(X), ncol = length(split_names),
                   dimnames = list(rownames(X), col_names))
  for (i in seq_along(split_names)) {
    newsig[, i] <- rowMeans(X[, split_names[[i]], drop = FALSE])
  }
  newsig
}

#' @keywords internal
expand_rows <- function(mat) {
  splits <- strsplit(rownames(mat), "\\|")
  idx <- rep(seq_len(nrow(mat)), times = lengths(splits))
  result <- mat[idx, , drop = FALSE]
  rownames(result) <- unlist(splits)
  result
}

#' @keywords internal
transferSymbol <- function(x) {
  alias2symbol <- read.csv(
    system.file("extdata", "NCBI_20251008_gene_result_alias2symbol.csv.gz", package = "SecAct"),
    as.is = TRUE
  )
  alias2symbol[is.na(alias2symbol[, "Alias"]), "Alias"] <- "NA"
  x[x %in% alias2symbol[, 1]] <- alias2symbol[
    match(x[x %in% alias2symbol[, 1]], alias2symbol[, 1]), 2
  ]
  x
}

#' @keywords internal
rm_duplicates <- function(mat) {
  gene_count <- table(rownames(mat))
  gene_dupl <- names(gene_count)[gene_count > 1]
  if (length(gene_dupl) > 0) {
    gene_unique <- names(gene_count)[gene_count == 1]
    gene_unique_index <- which(rownames(mat) %in% gene_unique)
    gene_dupl_index <- vapply(gene_dupl, function(gene) {
      idx <- which(rownames(mat) == gene)
      sums <- Matrix::rowSums(mat[idx, , drop = FALSE])
      idx[which.max(sums)]
    }, integer(1))
    mat <- mat[sort(c(gene_unique_index, gene_dupl_index)), ]
  }
  return(mat)
}

#' @keywords internal
sweep_sparse <- function(m, margin, stats, fun) {
  f <- match.fun(fun)
  if (margin == 1) {
    idx <- m@i + 1
  } else {
    if (inherits(m, "dgCMatrix")) {
      idx <- rep(1:m@Dim[2], diff(m@p))
    } else {
      idx <- m@j + 1
    }
  }
  m@x <- f(m@x, stats[idx])
  m
}

#' @keywords internal
.secact_preprocess <- function(Y, SigMat, is.group.sig, is.group.cor, rng_method) {
  if (is.matrix(SigMat) || is.data.frame(SigMat)) {
    X <- SigMat
  } else if (SigMat == "SecAct") {
    Xfile <- system.file("extdata", "SecAct.tsv.gz", package = "RidgeR")
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
  if (length(olp) < 2) stop("The overlapped genes between your expression matrix and the signature matrix are too few!")
  X <- as.matrix(X[olp, , drop = FALSE])
  Y <- as.matrix(Y[olp, , drop = FALSE])

  X <- scale(X); X[is.na(X)] <- 0
  Y <- scale(Y); Y[is.na(Y)] <- 0

  list(X = X, Y = Y, rng_int = rng_int)
}

#' @keywords internal
.secact_postprocess <- function(result, is.group.sig) {
  if (is.group.sig) {
    for (nm in c("beta", "se", "zscore", "pvalue")) {
      result[[nm]] <- expand_rows(result[[nm]])
    }
    ord <- sort(rownames(result$beta))
    for (nm in c("beta", "se", "zscore", "pvalue")) {
      result[[nm]] <- result[[nm]][ord, , drop = FALSE]
    }
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
# Platform-aware wrapper
# =========================================================

#' @title Secreted protein activity inference (Platform-aware wrapper)
#' @description Unified entry point that automatically selects the best backend
#'   for the current platform, or dispatches to a user-specified backend via the
#'   \code{method} parameter. Platform defaults:
#'   \itemize{
#'     \item \strong{Linux}: \code{SecAct.inference.Tcol.mt} — T-column permutation, multi-threaded (fastest)
#'     \item \strong{macOS}: \code{SecAct.inference.Yrow.st} — single-threaded (OpenMP unreliable without libomp, veclib)
#'     \item \strong{Windows}: \code{SecAct.inference.Tcol.mt} — T-column permutation, multi-threaded (OpenMP by Rtools)
#'   }
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: \code{"SecAct"} (bundled) or path to a tab-separated file.
#' @param lambda Ridge regularization parameter (default: 5e+05).
#' @param nrand Number of permutations (default: 1000).
#' @param ncores Number of CPU cores for multi-threaded variants (\code{NULL} = auto-detect).
#' @param rng_method RNG backend: \code{"srand"} (default, platform-dependent C stdlib)
#'   or \code{"gsl"} (cross-platform GSL MT19937, matches SecActPy).
#' @param method Backend to use. One of \code{"auto"} (default), \code{"Tcol.mt"},
#'   \code{"Tcol.st"}, \code{"Yrow.mt"}, \code{"Yrow.st"}, \code{"naive"}, or \code{"gsl.old"}.
#'   \code{"auto"} selects the platform default.
#' @param is.group.sig Logical; group correlated signatures before regression (default TRUE).
#' @param is.group.cor Correlation threshold for signature grouping (default 0.9).
#' @return A list with matrices \code{beta}, \code{se}, \code{zscore}, \code{pvalue}
#'   of dimension (proteins x samples).
#' @examples
#' \dontrun{
#' # Auto-select best backend for current platform
#' res <- SecAct.inference(expr.diff)
#'
#' # Force a specific backend
#' res <- SecAct.inference(expr.diff, method = "Tcol.mt", ncores = 8)
#' res <- SecAct.inference(expr.diff, method = "gsl.old")
#' }
#' @export
SecAct.inference <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                             ncores = NULL, rng_method = "srand", method = "auto",
                             is.group.sig = TRUE, is.group.cor = 0.9)
{
  method <- match.arg(method, c("auto", "Tcol.mt", "Tcol.st", "Yrow.mt", "Yrow.st", "naive", "gsl.old"))

  # Resolve "auto" to the platform-appropriate backend
  if (method == "auto") {
    is_mac     <- grepl("darwin",  R.version$os, ignore.case = TRUE)
    is_windows <- .Platform$OS.type == "windows"

    if (is_mac) {
      method <- "Yrow.st"   # macOS: OpenMP unreliable without libomp → single-threaded (w veclib blas)
      message("[SecAct.inference] Platform: macOS -> using Yrow.st")
    } else {
      method <- "Tcol.mt"   # Linux / Windows (Rtools): multi-threaded T-column permutation
      message("[SecAct.inference] Platform: ", ifelse(is_windows, "Windows", "Linux"), " -> using Tcol.mt")
    }
  }

  # Dispatch
  fn <- switch(method,
    Tcol.mt = SecAct.inference.Tcol.mt,
    Tcol.st = SecAct.inference.Tcol.st,
    Yrow.mt = SecAct.inference.Yrow.mt,
    Yrow.st = SecAct.inference.Yrow.st,
    naive   = SecAct.inference.naive,
    gsl.old = SecAct.inference.gsl.old
  )
  args <- list(Y = Y, SigMat = SigMat, lambda = lambda, nrand = nrand,
               rng_method = rng_method, is.group.sig = is.group.sig, is.group.cor = is.group.cor)
  if (method %in% c("Tcol.mt", "Yrow.mt")) args$ncores <- ncores
  do.call(fn, args)
}


# =========================================================
# Individual backends
# =========================================================

#' @title Secreted protein activity inference (Legacy)
#' @description Single-threaded legacy implementation using .C interface (32-bit indexing).
#'   Preserved for backward compatibility and macOS stability.
#' @param Y Gene expression matrix.
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @param rng_method RNG method: "srand" (default) or "gsl".
#' @param is.group.sig Logical; group correlated signatures before regression (default TRUE).
#' @param is.group.cor Correlation threshold for grouping (default 0.9).
#' @export
SecAct.inference.gsl.old <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     rng_method = "srand", is.group.sig = TRUE, is.group.cor = 0.9)
{
  pre <- .secact_preprocess(Y, SigMat, is.group.sig, is.group.cor, rng_method)

  n <- nrow(pre$X); p <- ncol(pre$X); m <- ncol(pre$Y)

  res <- .C("ridgeReg",
            X = as.double(t(pre$X)),
            Y = as.double(t(pre$Y)),
            as.integer(n), as.integer(p), as.integer(m),
            as.double(lambda), as.double(nrand),
            beta = double(p * m), se = double(p * m),
            zscore = double(p * m), pvalue = double(p * m),
            rng_method = as.integer(pre$rng_int)
  )

  result <- .secact_format_call(res, pre$X, pre$Y)
  .secact_postprocess(result, is.group.sig)
}

#' @title Secreted protein activity inference (Optimized dispatcher)
#' @description Dispatches to one of 4 optimized .Call variants. For platform-aware
#'   auto-selection use \code{SecAct.inference()}.
#' @param Y Gene expression matrix.
#' @param SigMat Secreted protein signature matrix.
#' @param lambda Penalty factor.
#' @param nrand Number of randomizations.
#' @param ncores Number of cores (NULL = auto-detect; multi-threaded variants only).
#' @param rng_method RNG method: "srand" (default) or "gsl".
#' @param method Variant: "Tcol.mt" (default), "Tcol.st", "Yrow.mt", or "Yrow.st".
#' @param is.group.sig Logical; group correlated signatures before regression (default TRUE).
#' @param is.group.cor Correlation threshold for grouping (default 0.9).
#' @export
SecAct.inference.gsl.new <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     ncores = NULL, rng_method = "srand", method = "Tcol.mt",
                                     is.group.sig = TRUE, is.group.cor = 0.9)
{
  method <- match.arg(method, c("Tcol.mt", "Tcol.st", "Yrow.mt", "Yrow.st"))
  SecAct.inference(Y, SigMat = SigMat, lambda = lambda, nrand = nrand,
                   ncores = ncores, rng_method = rng_method, method = method,
                   is.group.sig = is.group.sig, is.group.cor = is.group.cor)
}


#' @keywords internal
.secact_call_backend <- function(Y, SigMat, lambda, nrand, ncores, rng_method,
                                  is.group.sig, is.group.cor, c_symbol) {
  pre <- .secact_preprocess(Y, SigMat, is.group.sig, is.group.cor, rng_method)
  if (is.null(ncores)) ncores <- parallel::detectCores(logical = FALSE)
  res <- .Call(c_symbol, pre$X, pre$Y,
               as.numeric(lambda), as.integer(nrand),
               as.integer(ncores), as.integer(pre$rng_int))
  result <- .secact_format_call(res, pre$X, pre$Y)
  .secact_postprocess(result, is.group.sig)
}

# =========================================================
# 5 Optimized Variants
# =========================================================

#' @title Y-row permutation, single-threaded (.Call)
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param rng_method RNG method: "srand" (default) or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.Yrow.st <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     rng_method = "srand", is.group.sig = TRUE, is.group.cor = 0.9)
{
  SecAct.inference.Yrow.mt(Y, SigMat = SigMat, lambda = lambda, nrand = nrand,
                           ncores = 1L, rng_method = rng_method,
                           is.group.sig = is.group.sig, is.group.cor = is.group.cor)
}

#' @title Y-row permutation, multi-threaded (.Call)
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param ncores Number of CPU cores (NULL = auto-detect).
#' @param rng_method RNG method: "srand" (default) or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.Yrow.mt <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     ncores = NULL, rng_method = "srand",
                                     is.group.sig = TRUE, is.group.cor = 0.9)
{
  .secact_call_backend(Y, SigMat, lambda, nrand, ncores, rng_method,
                       is.group.sig, is.group.cor, "ridgeRegFast_interface")
}

#' @title T-column permutation, single-threaded (.Call)
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param rng_method RNG method: "srand" (default) or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.Tcol.st <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     rng_method = "srand", is.group.sig = TRUE, is.group.cor = 0.9)
{
  SecAct.inference.Tcol.mt(Y, SigMat = SigMat, lambda = lambda, nrand = nrand,
                           ncores = 1L, rng_method = rng_method,
                           is.group.sig = is.group.sig, is.group.cor = is.group.cor)
}

#' @title T-column permutation, multi-threaded (.Call)
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param ncores Number of CPU cores (NULL = auto-detect).
#' @param rng_method RNG method: "srand" (default) or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.Tcol.mt <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                     ncores = NULL, rng_method = "srand",
                                     is.group.sig = TRUE, is.group.cor = 0.9)
{
  .secact_call_backend(Y, SigMat, lambda, nrand, ncores, rng_method,
                       is.group.sig, is.group.cor, "ridgeRegFastTcol_interface")
}

#' @title Pure R naive implementation with C-generated permutation table
#' @param Y Gene expression matrix (genes x samples).
#' @param SigMat Signature matrix: "SecAct" (bundled) or path to file.
#' @param lambda Ridge regularization parameter.
#' @param nrand Number of permutations.
#' @param rng_method RNG method: "srand" (default) or "gsl".
#' @param is.group.sig Group correlated signatures before regression.
#' @param is.group.cor Correlation threshold for grouping.
#' @export
SecAct.inference.naive <- function(Y, SigMat = "SecAct", lambda = 5e+05, nrand = 1000,
                                   rng_method = "srand", is.group.sig = TRUE, is.group.cor = 0.9)
{
  pre <- .secact_preprocess(Y, SigMat, is.group.sig, is.group.cor, rng_method)
  X <- pre$X
  Y <- pre$Y

  n <- nrow(X); p <- ncol(X); m <- ncol(Y)

  # 1. Projection matrix: T = (X'X + lambda*I)^-1 X'
  T_mat <- solve(crossprod(X) + lambda * diag(p), t(X))  # (p x n)

  # 2. Observed beta
  beta <- T_mat %*% Y  # (p x m)

  # 3. Permutation table from C (same RNG sequence as C backends)
  perm_table <- .Call("generate_perm_table",
                      as.integer(n),
                      as.integer(nrand),
                      as.integer(pre$rng_int))
  # perm_table: (nrand x n) integer matrix, 0-indexed

  # 4. Permutation loop (T-column approach)
  aver          <- matrix(0.0, nrow = p, ncol = m)
  aver_sq       <- matrix(0.0, nrow = p, ncol = m)
  pvalue_counts <- matrix(0.0, nrow = p, ncol = m)
  abs_beta      <- abs(beta)

  for (i in seq_len(nrand)) {
    fwd_perm <- perm_table[i, ] + 1L
    inv_perm <- integer(n)
    inv_perm[fwd_perm] <- seq_len(n)
    beta_rand     <- T_mat[, inv_perm, drop = FALSE] %*% Y
    pvalue_counts <- pvalue_counts + (abs(beta_rand) >= abs_beta) * 1.0
    aver          <- aver    + beta_rand
    aver_sq       <- aver_sq + beta_rand^2
  }

  # 5. Finalize
  mean_rand <- aver / nrand
  var_rand  <- (aver_sq / nrand) - mean_rand^2
  var_rand[var_rand < 0] <- 0
  se     <- sqrt(var_rand)
  zscore <- ifelse(se > 1e-12, (beta - mean_rand) / se, 0.0)
  pvalue <- (pvalue_counts + 1.0) / (nrand + 1.0)

  result <- list(beta = beta, se = se, zscore = zscore, pvalue = pvalue)
  for (nm in names(result)) dimnames(result[[nm]]) <- list(colnames(X), colnames(Y))
  .secact_postprocess(result, is.group.sig)
}


# =========================================================
# SecAct activity inference wrappers
# =========================================================

#' @title Secreted protein activity inference
#' @description Infer the signaling activity of over 1000 secreted proteins from gene expression profiles.
#' @param inputProfile Gene expression matrix with gene symbol (row) x sample (column).
#' @param inputProfile_control Gene expression matrix with gene symbol (row) x sample (column).
#' @param is.differential A logical flag indicating whether inputProfile has been differential profiles against the control (Default: FALSE).
#' @param is.paired A logical flag indicating whether you want a paired operation of differential profiles between inputProfile and inputProfile_control if samples in inputProfile and inputProfile_control are paired (Default: FALSE).
#' @param is.singleSampleLevel A logical flag indicating whether to calculate activity change for each single sample between inputProfile and inputProfile_control (Default: FALSE). If FALSE, calculate the overall activity change between two phenotypes.
#' @param sigMatrix Secreted protein signature matrix. Could be "SecAct", "CytoSig", "SecAct-Breast", "SecAct-Colorectal", "SecAct-Glioblastoma", "SecAct-Kidney", "SecAct-Liver", "SecAct-Lung-Adeno", "SecAct-Ovarian", "SecAct-Pancreatic", "SecAct-Prostate". SecAct signatures were derived from all cancer ST samples; SecAct-XXX signatures were derived from XXX cancer ST samples.
#' @param is.filter.sig A logical flag indicating whether to filter the secreted protein signatures based on the genes from inputProfile (Default: FALSE). Because some sequencing platforms (e.g., CosMx) cover only a subset of secreted proteins, setting this option to TRUE restricts the activity inference on those proteins.
#' @param is.group.sig A logical flag indicating whether to group similar signatures (Default: TRUE). Many secreted proteins, such as cytokines with similar cell surface receptors and downstream pathways, have cellular effects that appear redundant within a cellular context. When enabled, this option clusters secreted proteins based on Pearson correlations among their composite signatures. The output still reports activity estimates for all secreted proteins prior to clustering. Secreted proteins assigned to the same non-redundant cluster share the same inferred activity.
#' @param is.group.cor A numeric value specifying the correlation cutoff used to define similar signatures (Default: 0.90). When r > 0.90, 1,170 secreted protein signatures are grouped into 657 non-redundant signature groups.
#' @param lambda Penalty factor in the ridge regression. If NULL, lambda will be assigned as 5e+05 or 10000 when sigMatrix = "SecAct" or "CytoSig", respectively.
#' @param nrand Number of randomization in the permutation test, with a default value of 1000.
#' @param ncores Number of CPU cores for multi-threaded variants (NULL = auto-detect).
#' @param rng_method RNG backend: "srand" (default) or "gsl".
#' @param method Backend to use: "auto" (default), "Tcol.mt", "Tcol.st", "Yrow.mt", "Yrow.st", "naive", or "gsl.old".
#' @return
#'
#' A list with four items, each is a matrix.
#' beta: regression coefficients
#' se: standard errors of coefficients
#' zscore: beta/se
#' pvalue: statistical significance
#'
#' @rdname SecAct.activity.inference
#' @export
#'
SecAct.activity.inference <- function(
    inputProfile,
    inputProfile_control = NULL,
    is.differential = FALSE,
    is.paired = FALSE,
    is.singleSampleLevel = FALSE,
    sigMatrix = "SecAct",
    is.filter.sig = FALSE,
    is.group.sig = TRUE,
    is.group.cor = 0.9,
    lambda = 5e+05,
    nrand = 1000,
    ncores = NULL,
    rng_method = "srand",
    method = "auto"
)
{
  if (inherits(inputProfile, "SpaCET")) {
    stop("Please use 'SecAct.activity.inference.ST'.")
  }
  if (inherits(inputProfile, "Seurat")) {
    stop("Please use 'SecAct.activity.inference.scRNAseq'.")
  }

  # --- Prepare Y (differential expression profiles) ---
  if (is.differential) {
    Y <- inputProfile
    if (ncol(Y) == 1) colnames(Y) <- "Change"
  } else {
    if (is.null(inputProfile_control)) {
      if (ncol(inputProfile) == 1) stop("Please include at least two samples in 'inputProfile' or set 'inputProfile_control'.")
      Y <- inputProfile - rowMeans(inputProfile)
    } else {
      if (is.paired) {
        Y <- inputProfile - inputProfile_control[, colnames(inputProfile), drop = FALSE]
      } else {
        Y <- inputProfile - rowMeans(inputProfile_control)
      }
      if (!is.singleSampleLevel) {
        Y <- matrix(rowMeans(Y), ncol = 1, dimnames = list(rownames(Y), "Change"))
      }
    }
  }

  # --- Load signature matrix ---
  if (sigMatrix == "SecAct") {
    Xfile <- system.file("extdata", "SecAct.tsv.gz", package = "RidgeR")
    X <- read.table(Xfile, sep = "\t", check.names = FALSE)
    if (is.null(lambda)) lambda <- 5e+05
  } else if (grepl("SecAct-", sigMatrix, fixed = TRUE)) {
    Xfile <- paste0("https://hpc.nih.gov/~Jiang_Lab/SecAct_Package/", sigMatrix, "_filterByPan_ds3_vst.tsv")
    X <- read.table(Xfile, sep = "\t", check.names = FALSE)
    if (is.null(lambda)) lambda <- 5e+05
  } else if (sigMatrix == "CytoSig") {
    Xfile <- "https://raw.githubusercontent.com/data2intelligence/CytoSig/refs/heads/master/CytoSig/signature.centroid"
    X <- read.table(Xfile, sep = "\t", check.names = FALSE)
    if (is.null(lambda)) lambda <- 10000
  } else {
    X <- read.table(sigMatrix, sep = "\t", check.names = FALSE)
  }

  # --- Filter signatures ---
  if (is.filter.sig) {
    X <- X[, colnames(X) %in% row.names(Y)]
  }

  # --- Call platform-aware inference ---
  SecAct.inference(
    Y = Y,
    SigMat = X,
    lambda = lambda,
    nrand = nrand,
    ncores = ncores,
    rng_method = rng_method,
    method = method,
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor
  )
}


#' @keywords internal
.preprocess_spacet_counts <- function(spacet_obj, scale.factor) {
  expr <- spacet_obj@input$counts
  expr <- expr[Matrix::rowSums(expr) > 0, ]
  rownames(expr) <- transferSymbol(rownames(expr))
  expr <- rm_duplicates(expr)
  stats <- Matrix::colSums(expr)
  expr <- sweep_sparse(expr, 2, stats, "/")
  expr@x <- expr@x * scale.factor
  expr@x <- log2(expr@x + 1)
  expr
}

#' @title Spot activity inference from spatial data
#' @description Calculate secreted protein signaling activity of spots from spatial transcriptomics data.
#' @param inputProfile A SpaCET object.
#' @param inputProfile_control A SpaCET object.
#' @param scale.factor Sets the scale factor for spot-level normalization.
#' @param sigMatrix Secreted protein signature matrix. Could be "SecAct", "CytoSig", "SecAct-Breast", "SecAct-Colorectal", "SecAct-Glioblastoma", "SecAct-Kidney", "SecAct-Liver", "SecAct-Lung-Adeno", "SecAct-Ovarian", "SecAct-Pancreatic", "SecAct-Prostate". SecAct signatures were derived from all cancer ST samples; SecAct-XXX signatures were derived from XXX cancer ST samples.
#' @param is.filter.sig A logical flag indicating whether to filter the secreted protein signatures based on the genes from inputProfile (Default: FALSE). Because some sequencing platforms (e.g., CosMx) cover only a subset of secreted proteins, setting this option to TRUE restricts the activity inference on those proteins.
#' @param is.group.sig A logical flag indicating whether to group similar signatures (Default: TRUE). Many secreted proteins, such as cytokines with similar cell surface receptors and downstream pathways, have cellular effects that appear redundant within a cellular context. When enabled, this option clusters secreted proteins based on Pearson correlations among their composite signatures. The output still reports activity estimates for all secreted proteins prior to clustering. Secreted proteins assigned to the same non-redundant cluster share the same inferred activity.
#' @param is.group.cor A numeric value specifying the correlation cutoff used to define similar signatures (Default: 0.90). When r > 0.90, 1,170 secreted protein signatures are grouped into 657 non-redundant signature groups.
#' @param lambda Penalty factor in the ridge regression. If NULL, lambda will be assigned as 5e+05 or 10000 when sigMatrix = "SecAct" or "CytoSig", respectively.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @param ncores Number of CPU cores for multi-threaded variants (NULL = auto-detect).
#' @param rng_method RNG backend: "srand" (default) or "gsl".
#' @param method Backend to use: "auto" (default), "Tcol.mt", "Tcol.st", "Yrow.mt", "Yrow.st", "naive", or "gsl.old".
#' @return A SpaCET object.
#' @rdname SecAct.activity.inference.ST
#' @export
#'
SecAct.activity.inference.ST <- function(
    inputProfile,
    inputProfile_control = NULL,
    scale.factor = 1e+05,
    sigMatrix = "SecAct",
    is.filter.sig = FALSE,
    is.group.sig = TRUE,
    is.group.cor = 0.9,
    lambda = 5e+05,
    nrand = 1000,
    ncores = NULL,
    rng_method = "srand",
    method = "auto"
)
{
  if (!inherits(inputProfile, "SpaCET")) {
    stop("Please input a SpaCET object.")
  }

  expr <- .preprocess_spacet_counts(inputProfile, scale.factor)

  if (is.null(inputProfile_control)) {
    expr.diff <- expr - Matrix::rowMeans(expr)
  } else {
    expr_control <- .preprocess_spacet_counts(inputProfile_control, scale.factor)
    olp <- intersect(rownames(expr), rownames(expr_control))
    expr.diff <- expr[olp, ] - Matrix::rowMeans(expr_control[olp, ])
  }

  res <- SecAct.activity.inference(
    inputProfile = expr.diff,
    is.differential = TRUE,
    sigMatrix = sigMatrix,
    is.filter.sig = is.filter.sig,
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor,
    lambda = lambda,
    nrand = nrand,
    ncores = ncores,
    rng_method = rng_method,
    method = method
  )

  inputProfile@results$SecAct_output$SecretedProteinActivity <- res

  inputProfile
}


#' @title Cell state activity inference from single cell data
#' @description Calculate secreted protein signaling activity of cell states from single cell RNA-Sequencing data.
#' @param inputProfile A Seurat object.
#' @param cellType_meta Column name in meta data that includes cell-type annotations.
#' @param is.singleCellLevel A logical flag indicating whether to calculate for each single cell (Default: FALSE).
#' @param sigMatrix Secreted protein signature matrix. Could be "SecAct", "CytoSig", "SecAct-Breast", "SecAct-Colorectal", "SecAct-Glioblastoma", "SecAct-Kidney", "SecAct-Liver", "SecAct-Lung-Adeno", "SecAct-Ovarian", "SecAct-Pancreatic", "SecAct-Prostate". SecAct signatures were derived from all cancer ST samples; SecAct-XXX signatures were derived from XXX cancer ST samples.
#' @param is.filter.sig A logical flag indicating whether to filter the secreted protein signatures based on the genes from inputProfile (Default: FALSE). Because some sequencing platforms (e.g., CosMx) cover only a subset of secreted proteins, setting this option to TRUE restricts the activity inference on those proteins.
#' @param is.group.sig A logical flag indicating whether to group similar signatures (Default: TRUE). Many secreted proteins, such as cytokines with similar cell surface receptors and downstream pathways, have cellular effects that appear redundant within a cellular context. When enabled, this option clusters secreted proteins based on Pearson correlations among their composite signatures. The output still reports activity estimates for all secreted proteins prior to clustering. Secreted proteins assigned to the same non-redundant cluster share the same inferred activity.
#' @param is.group.cor A numeric value specifying the correlation cutoff used to define similar signatures (Default: 0.90). When r > 0.90, 1,170 secreted protein signatures are grouped into 657 non-redundant signature groups.
#' @param lambda Penalty factor in the ridge regression. If NULL, lambda will be assigned as 5e+05 or 10000 when sigMatrix = "SecAct" or "CytoSig", respectively.
#' @param nrand Number of randomization in the permutation test, with a default value 1000.
#' @param ncores Number of CPU cores for multi-threaded variants (NULL = auto-detect).
#' @param rng_method RNG backend: "srand" (default) or "gsl".
#' @param method Backend to use: "auto" (default), "Tcol.mt", "Tcol.st", "Yrow.mt", "Yrow.st", "naive", or "gsl.old".
#' @return A Seurat object.
#' @rdname SecAct.activity.inference.scRNAseq
#' @export
#'
SecAct.activity.inference.scRNAseq <- function(
    inputProfile,
    cellType_meta,
    is.singleCellLevel = FALSE,
    sigMatrix = "SecAct",
    is.filter.sig = FALSE,
    is.group.sig = TRUE,
    is.group.cor = 0.9,
    lambda = 5e+05,
    nrand = 1000,
    ncores = NULL,
    rng_method = "srand",
    method = "auto"
)
{
  if (!inherits(inputProfile, "Seurat")) {
    stop("Please input a Seurat object.")
  }

  if (inherits(inputProfile@assays$RNA, "Assay5")) {
    counts <- inputProfile@assays$RNA@layers$counts
    colnames(counts) <- rownames(inputProfile@assays$RNA@cells)
    rownames(counts) <- rownames(inputProfile@assays$RNA@features)
  } else {
    counts <- inputProfile@assays$RNA@counts
  }

  rownames(counts) <- transferSymbol(rownames(counts))
  counts <- rm_duplicates(counts)

  if (!is.singleCellLevel) {
    cellType_vec <- inputProfile@meta.data[, cellType_meta]

    # generate pseudo bulk
    cell_types <- sort(unique(cellType_vec))
    expr <- matrix(NA_real_, nrow = nrow(counts), ncol = length(cell_types),
                   dimnames = list(rownames(counts), cell_types))
    for (i in seq_along(cell_types)) {
      expr[, i] <- Matrix::rowSums(counts[, cellType_vec == cell_types[i], drop = FALSE])
    }

    # normalize to TPM
    expr <- sweep(expr, 2, colSums(expr), "/") * 1e6
  } else {
    expr <- sweep(counts, 2, Matrix::colSums(counts), "/") * 1e5
  }
  rm(counts); gc()

  # transform to log space
  expr <- log2(expr + 1)

  # normalized with the control samples
  expr.diff <- expr - Matrix::rowMeans(expr)

  rm(expr); gc()

  inputProfile@misc$SecAct_output$SecretedProteinActivity <-
    SecAct.activity.inference(
      inputProfile = expr.diff,
      is.differential = TRUE,
      sigMatrix = sigMatrix,
      is.filter.sig = is.filter.sig,
      is.group.sig = is.group.sig,
      is.group.cor = is.group.cor,
      lambda = lambda,
      nrand = nrand,
      ncores = ncores,
      rng_method = rng_method,
      method = method
    )

  inputProfile
}
