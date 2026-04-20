# Canonical ridge() + ridge_batch() entry points.
#
# These exist to match the shared accelerator API used by
# RidgeFast::ridge / RidgeCuda::ridge, so ridge-bench's unified
# dispatcher and the cross-package bit-parity test can call RidgeR
# through the same signature as every other end-product package.
#
# They are *additions*, not replacements: SecAct.inference.*
# variants (Yrow/Tcol × st/mt, naive, gsl.old, gsl.new) remain the
# research-facing API and are unchanged.

#' Canonical ridge regression with permutation testing
#'
#' Direct dispatch to RidgeR's T-column C core. Matches
#' \code{RidgeFast::ridge} byte-for-byte in its permutation stream
#' when \code{rng_method = "srand"} + same \code{seed}.
#'
#' @param X Numeric matrix, n x p, column-scaled signature matrix.
#'   User is responsible for pre-scaling (column z-score). No
#'   internal scaling, no SigMat loading, no signature grouping.
#' @param Y Numeric matrix, n x m, column-scaled expression matrix.
#' @param lambda Ridge penalty (default 5e+05).
#' @param nrand Number of permutations (default 1000).
#' @param ncores OpenMP threads (default 1). Use 1 for
#'   bit-reproducible output.
#' @param rng_method \code{"srand"} (default; matches original SecAct
#'   C behavior) or \code{"gsl"} (GSL MT19937, cross-platform
#'   reproducible).
#' @param seed Integer seed for the RNG (default 0).
#' @return A list with four p-by-m matrices (\code{beta}, \code{se},
#'   \code{zscore}, \code{pvalue}).
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200; p <- 20; m <- 10
#' X <- scale(matrix(rnorm(n * p), n, p)); colnames(X) <- paste0("sig", 1:p)
#' Y <- scale(matrix(rnorm(n * m), n, m)); colnames(Y) <- paste0("s", 1:m)
#' res <- ridge(X, Y, lambda = 1, nrand = 100)
#' str(res)
#' }
#' @useDynLib RidgeR, .registration = TRUE
#' @export
ridge <- function(X, Y, lambda = 5e+05, nrand = 1000L,
                  ncores = 1L, rng_method = "srand", seed = 0L) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.numeric(X) || !is.numeric(Y)) stop("X and Y must be numeric.")
  if (nrow(X) != nrow(Y)) stop("nrow(X) must equal nrow(Y).")
  if (!is.finite(lambda) || lambda < 0) stop("lambda must be non-negative finite.")
  if (nrand < 1L) stop("nrand must be >= 1.")
  if (ncores < 1L) stop("ncores must be >= 1.")

  rng_method <- match.arg(rng_method, c("srand", "gsl"))
  rng_int <- if (rng_method == "srand") 0L else 1L

  storage.mode(X) <- "double"
  storage.mode(Y) <- "double"

  res <- .Call("ridgeRegFastTcol_interface",
               X, Y,
               as.double(lambda), as.integer(nrand),
               as.integer(ncores), rng_int, as.integer(seed),
               PACKAGE = "RidgeR")

  p <- ncol(X); m <- ncol(Y)
  dn <- list(colnames(X), colnames(Y))
  list(
    beta   = matrix(res$beta,   p, m, byrow = TRUE, dimnames = dn),
    se     = matrix(res$se,     p, m, byrow = TRUE, dimnames = dn),
    zscore = matrix(res$zscore, p, m, byrow = TRUE, dimnames = dn),
    pvalue = matrix(res$pvalue, p, m, byrow = TRUE, dimnames = dn)
  )
}


#' Batched canonical ridge regression
#'
#' Memory-efficient wrapper around \code{\link{ridge}} that processes
#' \code{Y} in column-batches. Supports in-memory matrices, HDF5
#' files, and user-supplied readers; results can be accumulated in
#' memory or streamed to HDF5. API mirrors
#' \code{RidgeFast::ridge_batch}.
#'
#' Observed \eqn{\beta} and permutation statistics are independent
#' across columns of \code{Y}, so batched output is identical (up to
#' FP accumulation order) to a single \code{ridge()} call at
#' \code{ncores = 1}. Each batch rebuilds the projection + permutation
#' table — that overhead is small for typical p.
#'
#' @param X Numeric matrix, n x p, column-scaled signature matrix.
#' @param Y Either (a) a dense matrix n x m, (b) a path to an HDF5
#'   file containing a dataset named \code{"Y"} of shape n x m
#'   (requires the \pkg{rhdf5} package), or (c) \code{NULL} together
#'   with a \code{reader} callback.
#' @param lambda Ridge penalty.
#' @param nrand Number of permutations.
#' @param ncores OpenMP threads per batch.
#' @param rng_method \code{"srand"} or \code{"gsl"}.
#' @param seed Integer seed.
#' @param batch_size Number of Y-columns per batch.
#' @param reader Optional \code{function(start, end)} returning the
#'   n x (end - start + 1) slice of \code{Y}.
#' @param n_samples Required when \code{reader} is used.
#' @param output_h5 Optional path to an HDF5 file; if supplied the
#'   four p x m result matrices stream to datasets
#'   \code{"beta" / "se" / "zscore" / "pvalue"}.
#' @param verbose Print per-batch progress.
#' @return Either a list of four p x m matrices or, when
#'   \code{output_h5} is set, metadata about the written file.
#' @seealso \code{\link{ridge}}
#' @export
ridge_batch <- function(X, Y, lambda = 5e+05, nrand = 1000L,
                        ncores = 1L, rng_method = "srand", seed = 0L,
                        batch_size = 5000L,
                        reader = NULL, n_samples = NULL,
                        output_h5 = NULL, verbose = FALSE) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  storage.mode(X) <- "double"
  n <- nrow(X); p <- ncol(X)
  sig_names <- colnames(X)
  if (is.null(sig_names)) sig_names <- paste0("sig", seq_len(p))

  batch_size <- as.integer(batch_size)
  if (is.na(batch_size) || batch_size < 1L) stop("batch_size must be >= 1.")

  src <- .resolve_y_source_canonical(Y, reader, n_samples, n)
  m <- src$m
  samp_names <- src$samp_names
  reader_fn <- src$reader_fn

  h5_out <- !is.null(output_h5)
  if (h5_out) {
    if (!requireNamespace("rhdf5", quietly = TRUE))
      stop("'output_h5' requires the 'rhdf5' package.")
    if (file.exists(output_h5)) file.remove(output_h5)
    rhdf5::h5createFile(output_h5)
    chunk_m <- min(batch_size, m)
    for (nm in c("beta", "se", "zscore", "pvalue")) {
      rhdf5::h5createDataset(output_h5, nm, dims = c(p, m),
                             storage.mode = "double",
                             chunk = c(p, chunk_m))
    }
    rhdf5::h5write(sig_names, output_h5, "signature_names")
    rhdf5::h5write(samp_names, output_h5, "sample_names")
  } else {
    out_beta   <- matrix(0, p, m, dimnames = list(sig_names, samp_names))
    out_se     <- matrix(0, p, m, dimnames = list(sig_names, samp_names))
    out_zscore <- matrix(0, p, m, dimnames = list(sig_names, samp_names))
    out_pvalue <- matrix(0, p, m, dimnames = list(sig_names, samp_names))
  }

  num_batches <- as.integer(ceiling(m / batch_size))
  for (b in seq_len(num_batches)) {
    s <- (b - 1L) * batch_size + 1L
    e <- min(as.integer(b * batch_size), m)
    if (verbose) message(sprintf("[RidgeR::ridge_batch] batch %d/%d (cols %d-%d)",
                                 b, num_batches, s, e))

    Y_batch <- reader_fn(s, e)
    if (!is.matrix(Y_batch)) Y_batch <- as.matrix(Y_batch)
    storage.mode(Y_batch) <- "double"
    if (nrow(Y_batch) != n)
      stop(sprintf("reader returned %d rows, expected %d.", nrow(Y_batch), n))
    if (ncol(Y_batch) != (e - s + 1L))
      stop(sprintf("reader returned %d cols for batch (%d-%d), expected %d.",
                   ncol(Y_batch), s, e, e - s + 1L))
    if (is.null(colnames(Y_batch))) colnames(Y_batch) <- samp_names[s:e]

    res <- ridge(X, Y_batch, lambda = lambda, nrand = nrand,
                 ncores = ncores, rng_method = rng_method, seed = seed)

    if (h5_out) {
      for (nm in c("beta", "se", "zscore", "pvalue")) {
        rhdf5::h5write(res[[nm]], output_h5, nm, index = list(NULL, s:e))
      }
    } else {
      out_beta[, s:e]   <- res$beta
      out_se[, s:e]     <- res$se
      out_zscore[, s:e] <- res$zscore
      out_pvalue[, s:e] <- res$pvalue
    }
  }

  if (h5_out) {
    rhdf5::h5closeAll()
    invisible(list(path = output_h5, p = p, m = m, num_batches = num_batches))
  } else {
    list(beta = out_beta, se = out_se, zscore = out_zscore, pvalue = out_pvalue)
  }
}

.resolve_y_source_canonical <- function(Y, reader, n_samples, n) {
  if (is.character(Y) && length(Y) == 1L) {
    if (!requireNamespace("rhdf5", quietly = TRUE))
      stop("Reading Y from an HDF5 file requires the 'rhdf5' package.")
    info <- rhdf5::h5ls(Y)
    y_row <- info[info$name == "Y", ]
    if (nrow(y_row) == 0L) stop("HDF5 file '", Y, "' must contain a dataset named 'Y'.")
    dims <- as.integer(strsplit(y_row$dim[1], " x ", fixed = TRUE)[[1]])
    if (length(dims) != 2L) stop("'Y' dataset must be 2-dimensional.")
    if (dims[1] != n) stop(sprintf("HDF5 Y rows (%d) do not match nrow(X) (%d).", dims[1], n))
    m <- dims[2]
    samp_names <- if ("sample_names" %in% info$name) {
      as.character(rhdf5::h5read(Y, "sample_names"))
    } else paste0("s", seq_len(m))
    path <- Y
    list(m = m, samp_names = samp_names,
         reader_fn = function(s, e) {
           rhdf5::h5read(path, "Y", index = list(NULL, s:e))
         })
  } else if (!is.null(reader)) {
    if (is.null(n_samples)) stop("'reader' requires 'n_samples'.")
    m <- as.integer(n_samples)
    list(m = m, samp_names = paste0("s", seq_len(m)), reader_fn = reader)
  } else {
    if (!is.matrix(Y)) Y <- as.matrix(Y)
    if (nrow(Y) != n) stop(sprintf("nrow(Y) (%d) != nrow(X) (%d).", nrow(Y), n))
    m <- ncol(Y)
    samp_names <- if (is.null(colnames(Y))) paste0("s", seq_len(m)) else colnames(Y)
    list(m = m, samp_names = samp_names,
         reader_fn = function(s, e) Y[, s:e, drop = FALSE])
  }
}
