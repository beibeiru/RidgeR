#' =============================================================================
#' RidgeR: Secreted Protein Activity Inference
#' =============================================================================
#'
#' This module provides R wrappers for ridge regression-based activity inference.
#' Multiple implementations are available with different performance characteristics:
#'
#'   - legacy: Original .C interface (single-threaded, Y-permutation)
#'   - styp:   .Call interface, Single-Thread Y-Permutation
#'   - sttp:   .Call interface, Single-Thread T-Permutation (cache-friendly)
#'   - mtyp:   .Call interface, Multi-Thread Y-Permutation (OpenMP)
#'   - mttp:   .Call interface, Multi-Thread T-Permutation (OpenMP, cache-friendly)
#'   - r:      Pure R implementation (no GSL required)
#'
#' @name RidgeR-package
#' @docType package
NULL


# ==============================================================================
# SECTION 1: INTERNAL HELPER FUNCTIONS
# ==============================================================================

#' Load signature matrix from file or package data
#'
#' @param SigMat Either "SecAct" for default or path to custom signature matrix
#' @return Data frame containing the signature matrix
#' @keywords internal
.load_signature <- function(SigMat) {
    if (SigMat == "SecAct") {
        sig_path <- file.path(
            system.file(package = "RidgeR"),
            "extdata/SecAct.tsv.gz"
        )
        if (!file.exists(sig_path)) {
            stop("Default signature matrix not found at: ", sig_path)
        }
        read.table(sig_path, sep = "\t", check.names = FALSE)
    } else {
        if (!file.exists(SigMat)) {
            stop("Signature matrix file not found: ", SigMat)
        }
        read.table(SigMat, sep = "\t", check.names = FALSE)
    }
}


#' Expand rows with pipe-delimited names
#'
#' Expands rows where rownames contain "|" (grouped signatures) into
#' separate rows with duplicated values.
#'
#' @param mat Matrix with potentially grouped row names (e.g., "A|B|C")
#' @return Matrix with expanded rows
#' @keywords internal
.expand_rows <- function(mat) {
    row_names <- rownames(mat)
    
    # Find rows that need expansion (contain "|")
    needs_expansion <- grepl("\\|", row_names)
    
    if (!any(needs_expansion)) {
        return(mat)
    }
    
    # Build new matrix
    new_rows <- list()
    new_names <- character(0)
    
    for (i in seq_len(nrow(mat))) {
        if (needs_expansion[i]) {
            # Split the grouped name
            split_names <- strsplit(row_names[i], "\\|")[[1]]
            for (name in split_names) {
                new_rows[[length(new_rows) + 1]] <- mat[i, , drop = FALSE]
                new_names <- c(new_names, name)
            }
        } else {
            new_rows[[length(new_rows) + 1]] <- mat[i, , drop = FALSE]
            new_names <- c(new_names, row_names[i])
        }
    }
    
    # Combine into matrix
    result <- do.call(rbind, new_rows)
    rownames(result) <- new_names
    
    result
}


#' Group similar signatures by correlation
#'
#' Clusters signatures by Pearson correlation and groups those with
#' correlation above the threshold.
#'
#' @param X Signature matrix (genes x proteins)
#' @param cor_threshold Correlation threshold for grouping (default: 0.9)
#' @return Matrix with grouped signatures (averaged within groups)
#' @keywords internal
.group_signatures <- function(X, cor_threshold = 0.9) {
    # Calculate distance based on correlation
    dis <- as.dist(1 - cor(X, method = "pearson"))
    hc <- hclust(dis, method = "complete")
    
    # Cut tree at correlation threshold
    group_labels <- cutree(hc, h = 1 - cor_threshold)
    
    # Create new signature matrix with grouped signatures
    newsig <- data.frame(row.names = rownames(X))
    
    for (j in unique(group_labels)) {
        gene_groups <- names(group_labels)[group_labels == j]
        group_name <- paste0(gene_groups, collapse = "|")
        newsig[[group_name]] <- rowMeans(X[, gene_groups, drop = FALSE])
    }
    
    as.matrix(newsig)
}


#' Prepare and align matrices for ridge regression
#'
#' Aligns X and Y by common row names, scales both matrices,
#' and replaces NA values with zeros.
#'
#' @param X Signature matrix (genes x proteins)
#' @param Y Expression matrix (genes x samples)
#' @param sort_genes Logical; if TRUE, sort common genes alphabetically
#' @return List with aligned and scaled X and Y matrices
#' @keywords internal
.prepare_matrices <- function(X, Y, sort_genes = FALSE) {
    # Find overlapping genes
    common_genes <- intersect(rownames(Y), rownames(X))

    if (length(common_genes) == 0) {
        stop("No common genes found between signature and expression matrices")
    }

    # Optionally sort genes alphabetically
    if (sort_genes) {
        common_genes <- sort(common_genes)
    }

    # Subset and scale
    X <- scale(as.matrix(X[common_genes, , drop = FALSE]))
    Y <- scale(as.matrix(Y[common_genes, , drop = FALSE]))

    # Replace NAs with zeros (from scaling constant columns)
    X[is.na(X)] <- 0
    Y[is.na(Y)] <- 0

    list(X = X, Y = Y)
}


#' Format result vector as named matrix
#'
#' GSL stores matrices in row-major order, so we need byrow=TRUE
#' to correctly interpret the flat vector in R's column-major convention.
#'
#' @param v Flat result vector (from GSL row-major matrix)
#' @param X Signature matrix (for row names = proteins)
#' @param Y Expression matrix (for column names = samples)
#' @return Matrix with appropriate dimensions and names
#' @keywords internal
.format_result <- function(v, X, Y) {
    matrix(
        v,
        nrow = ncol(X),
        ncol = ncol(Y),
        byrow = TRUE,  # GSL uses row-major order
        dimnames = list(colnames(X), colnames(Y))
    )
}


#' Determine default number of cores
#'
#' @param ncores User-specified core count or NULL
#' @return Number of cores to use (at least 1)
#' @keywords internal
.default_ncores <- function(ncores) {
    if (is.null(ncores)) {
        max(1, parallel::detectCores(logical = FALSE) - 1)
    } else {
        as.integer(ncores)
    }
}


# ==============================================================================
# SECTION 2: LEGACY IMPLEMENTATION (.C INTERFACE)
# ==============================================================================

#' Secreted Protein Activity Inference (Legacy .C Version)
#'
#' Original implementation using the legacy .C interface.
#' This version is maintained for backward compatibility.
#'
#' @param Y Gene expression matrix (genes x samples)
#' @param SigMat Secreted protein signature matrix path or "SecAct" for default
#' @param lambda Ridge penalty factor (default: 5e+05)
#' @param nrand Number of permutations for significance testing (default: 1000)
#'
#' @return Named list with four matrices (proteins x samples):
#'   \itemize{
#'     \item \code{beta}: Ridge regression coefficients
#'     \item \code{se}: Standard errors from permutation null distribution
#'     \item \code{zscore}: Z-scores (beta normalized by permutation distribution)
#'     \item \code{pvalue}: Permutation-based p-values
#'   }
#'
#' @export
SecAct.inference.gsl.legacy <- function(Y,
                                         SigMat = "SecAct",
                                         lambda = 5e+05,
                                         nrand = 1000) {
    # Load and prepare data
    X <- .load_signature(SigMat)
    prep <- .prepare_matrices(X, Y)

    X <- prep$X
    Y <- prep$Y

    # Dimensions
    n <- as.integer(nrow(X))  # Number of genes
    p <- as.integer(ncol(X))  # Number of proteins
    m <- as.integer(ncol(Y))  # Number of samples
    len <- p * m

    # Call C function - NOTE: must transpose matrices for correct memory layout
    # R stores column-major, GSL expects row-major, so t(X) in R = X in GSL
    res <- .C(
        "ridgeReg",
        as.double(t(X)),   # Transpose for correct GSL interpretation
        as.double(t(Y)),   # Transpose for correct GSL interpretation
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

    # Format and return results - note byrow=TRUE to match original
    list(
        beta   = matrix(res$beta, byrow = TRUE, ncol = m,
                       dimnames = list(colnames(X), colnames(Y))),
        se     = matrix(res$se, byrow = TRUE, ncol = m,
                       dimnames = list(colnames(X), colnames(Y))),
        zscore = matrix(res$zscore, byrow = TRUE, ncol = m,
                       dimnames = list(colnames(X), colnames(Y))),
        pvalue = matrix(res$pvalue, byrow = TRUE, ncol = m,
                       dimnames = list(colnames(X), colnames(Y)))
    )
}


# ==============================================================================
# SECTION 3: SINGLE-THREADED IMPLEMENTATIONS (.Call INTERFACE)
# ==============================================================================

#' Secreted Protein Activity Inference (Single-threaded, Y-Permutation)
#'
#' Single-threaded implementation using the .Call interface with
#' row permutation of the response matrix Y.
#'
#' @inheritParams SecAct.inference.gsl.legacy
#'
#' @return Named list with beta, se, zscore, pvalue matrices (proteins x samples)
#'
#' @export
SecAct.inference.gsl.styp <- function(Y,
                                      SigMat = "SecAct",
                                      lambda = 5e+05,
                                      nrand = 1000) {
    # Load and prepare data
    X <- .load_signature(SigMat)
    prep <- .prepare_matrices(X, Y)

    # Call C function
    res <- .Call(
        "ridgeReg_styp_interface",
        prep$X,
        prep$Y,
        as.numeric(lambda),
        as.integer(nrand),
        PACKAGE = "RidgeR"
    )

    # Format and return results
    list(
        beta   = .format_result(res$beta,   prep$X, prep$Y),
        se     = .format_result(res$se,     prep$X, prep$Y),
        zscore = .format_result(res$zscore, prep$X, prep$Y),
        pvalue = .format_result(res$pvalue, prep$X, prep$Y)
    )
}


#' Secreted Protein Activity Inference (Single-threaded, T-Permutation)
#'
#' Single-threaded implementation using cache-friendly T column permutation.
#' This approach can be faster than Y-permutation for certain matrix sizes
#' due to better memory access patterns.
#'
#' @inheritParams SecAct.inference.gsl.legacy
#'
#' @return Named list with beta, se, zscore, pvalue matrices (proteins x samples)
#'
#' @export
SecAct.inference.gsl.sttp <- function(Y,
                                       SigMat = "SecAct",
                                       lambda = 5e+05,
                                       nrand = 1000) {
    # Load and prepare data
    X <- .load_signature(SigMat)
    prep <- .prepare_matrices(X, Y)

    # Call C function
    res <- .Call(
        "ridgeReg_sttp_interface",
        prep$X,
        prep$Y,
        as.numeric(lambda),
        as.integer(nrand),
        PACKAGE = "RidgeR"
    )

    # Format and return results
    list(
        beta   = .format_result(res$beta,   prep$X, prep$Y),
        se     = .format_result(res$se,     prep$X, prep$Y),
        zscore = .format_result(res$zscore, prep$X, prep$Y),
        pvalue = .format_result(res$pvalue, prep$X, prep$Y)
    )
}


# ==============================================================================
# SECTION 4: MULTI-THREADED IMPLEMENTATIONS (OpenMP)
# ==============================================================================

#' Secreted Protein Activity Inference (Multi-threaded, Y-Permutation)
#'
#' Fast multi-threaded implementation using OpenMP parallelization
#' with row permutation of the response matrix Y.
#'
#' @inheritParams SecAct.inference.gsl.legacy
#' @param ncores Number of CPU cores to use. Default (NULL) uses all
#'   available physical cores minus one.
#'
#' @return Named list with beta, se, zscore, pvalue matrices (proteins x samples)
#'
#' @export
SecAct.inference.gsl.mtyp <- function(Y,
                                      SigMat = "SecAct",
                                      lambda = 5e+05,
                                      nrand = 1000,
                                      ncores = NULL) {
    # Load and prepare data
    X <- .load_signature(SigMat)
    prep <- .prepare_matrices(X, Y)
    ncores <- .default_ncores(ncores)

    # Call C function
    res <- .Call(
        "ridgeReg_mtyp_interface",
        prep$X,
        prep$Y,
        as.numeric(lambda),
        as.integer(nrand),
        as.integer(ncores),
        PACKAGE = "RidgeR"
    )

    # Format and return results
    list(
        beta   = .format_result(res$beta,   prep$X, prep$Y),
        se     = .format_result(res$se,     prep$X, prep$Y),
        zscore = .format_result(res$zscore, prep$X, prep$Y),
        pvalue = .format_result(res$pvalue, prep$X, prep$Y)
    )
}


#' Secreted Protein Activity Inference (Multi-threaded, T-Permutation)
#'
#' Fast multi-threaded implementation using OpenMP parallelization
#' with cache-friendly T column permutation (scatter approach).
#'
#' @inheritParams SecAct.inference.gsl.mtyp
#'
#' @return Named list with beta, se, zscore, pvalue matrices (proteins x samples)
#'
#' @export
SecAct.inference.gsl.mttp <- function(Y,
                                       SigMat = "SecAct",
                                       lambda = 5e+05,
                                       nrand = 1000,
                                       ncores = NULL) {
    # Load and prepare data
    X <- .load_signature(SigMat)
    prep <- .prepare_matrices(X, Y)
    ncores <- .default_ncores(ncores)

    # Call C function
    res <- .Call(
        "ridgeReg_mttp_interface",
        prep$X,
        prep$Y,
        as.numeric(lambda),
        as.integer(nrand),
        as.integer(ncores),
        PACKAGE = "RidgeR"
    )

    # Format and return results
    list(
        beta   = .format_result(res$beta,   prep$X, prep$Y),
        se     = .format_result(res$se,     prep$X, prep$Y),
        zscore = .format_result(res$zscore, prep$X, prep$Y),
        pvalue = .format_result(res$pvalue, prep$X, prep$Y)
    )
}


# ==============================================================================
# SECTION 5: PURE R IMPLEMENTATION
# ==============================================================================

#' Secreted Protein Activity Inference (Pure R Implementation)
#'
#' Pure R implementation of ridge regression with permutation testing.
#' This version does not require GSL and is useful for debugging,
#' validation, and environments where GSL is unavailable.
#'
#' @inheritParams SecAct.inference.gsl.legacy
#' @param is.group.sig Logical indicating whether to group similar signatures
#'   (default: FALSE for this low-level function).
#' @param is.group.cor Correlation cutoff for signature grouping (default: 0.9).
#' @param sort_genes Logical; if TRUE, sort common genes alphabetically for
#'   reproducible ordering across different platforms (default: FALSE).
#'
#' @return Named list with beta, se, zscore, pvalue matrices (proteins x samples)
#'
#' @details
#' This implementation uses R's native linear algebra functions:
#' \itemize{
#'   \item \code{crossprod()} for X'X computation
#'   \item \code{chol()} for Cholesky decomposition
#'   \item \code{backsolve()} and \code{forwardsolve()} for solving linear systems
#' }
#'
#' Each permutation uses \code{set.seed(i)} where i is the permutation index,
#' ensuring reproducibility across runs.
#'
#' @examples
#' \dontrun{
#' # Create test data
#' set.seed(123)
#' Y <- matrix(rnorm(1000 * 10), nrow = 1000, ncol = 10)
#' rownames(Y) <- paste0("Gene", 1:1000)
#' colnames(Y) <- paste0("Sample", 1:10)
#'
#' # Run inference
#' res <- SecAct.inference.r(Y, nrand = 100)
#' head(res$zscore)
#' }
#'
#' @export
SecAct.inference.r <- function(Y,
                                SigMat = "SecAct",
                                lambda = 5e+05,
                                nrand = 1000,
                                is.group.sig = FALSE,
                                is.group.cor = 0.9,
                                sort_genes = FALSE) {
    # Load signature matrix
    X <- .load_signature(SigMat)
    
    # Group similar signatures if requested
    if (is.group.sig) {
        X <- .group_signatures(X, is.group.cor)
    }

    # Find overlapping genes
    olp <- intersect(rownames(Y), rownames(X))
    if (length(olp) < 2) {
        stop("Too few overlapping genes between expression and signature matrices!")
    }

    # Optionally sort genes alphabetically
    if (sort_genes) {
        olp <- sort(olp)
    }

    # Subset and scale
    X <- as.matrix(X[olp, , drop = FALSE])
    Y <- as.matrix(Y[olp, , drop = FALSE])
    X <- scale(X)
    Y <- scale(Y)

    # Handle NAs from scaling
    X[is.na(X)] <- 0
    Y[is.na(Y)] <- 0

    # Dimensions
    n <- nrow(Y)
    p <- ncol(X)
    m <- ncol(Y)

    # Compute ridge regression: beta = (X'X + lambda*I)^-1 * X' * Y
    # Use Cholesky decomposition for numerical stability
    A <- crossprod(X) + lambda * diag(p)  # X'X + lambda*I (SPD)
    R <- chol(A)                           # A = R'R

    # Solve for beta: R'R * beta = X'Y
    beta <- backsolve(R, forwardsolve(t(R), crossprod(X, Y)))

    # Permutation testing
    aver <- NULL
    aver_sq <- NULL
    pvalue <- NULL

    for (i in seq_len(nrand)) {
        set.seed(i)
        beta_rand <- backsolve(R, forwardsolve(t(R), crossprod(X, Y[sample.int(n), , drop = FALSE])))

        if (i == 1) {
            aver <- beta_rand
            aver_sq <- beta_rand^2
            pvalue <- (abs(beta_rand) >= abs(beta)) * 1.0
        } else {
            aver <- aver + beta_rand
            aver_sq <- aver_sq + beta_rand^2
            pvalue <- pvalue + (abs(beta_rand) >= abs(beta))
        }
    }

    # Finalize statistics
    aver <- aver / nrand
    aver_sq <- aver_sq / nrand
    se <- sqrt(aver_sq - aver * aver)
    zscore <- (beta - aver) / se
    zscore[!is.finite(zscore)] <- 0
    pvalue <- (pvalue + 1) / (nrand + 1)

    # Set dimension names
    rownames(beta) <- colnames(X)
    colnames(beta) <- colnames(Y)
    rownames(se) <- colnames(X)
    colnames(se) <- colnames(Y)
    rownames(zscore) <- colnames(X)
    colnames(zscore) <- colnames(Y)
    rownames(pvalue) <- colnames(X)
    colnames(pvalue) <- colnames(Y)

    # Expand grouped signatures back to individual rows
    if (is.group.sig) {
        beta <- .expand_rows(beta)
        se <- .expand_rows(se)
        zscore <- .expand_rows(zscore)
        pvalue <- .expand_rows(pvalue)

        # Sort by row name
        row_order <- sort(rownames(beta))
        beta <- beta[row_order, , drop = FALSE]
        se <- se[row_order, , drop = FALSE]
        zscore <- zscore[row_order, , drop = FALSE]
        pvalue <- pvalue[row_order, , drop = FALSE]
    }

    list(beta = beta, se = se, zscore = zscore, pvalue = pvalue)
}


# ==============================================================================
# SECTION 6: MAIN INFERENCE FUNCTION
# ==============================================================================

#' Secreted Protein Activity Inference (Full Pipeline)
#'
#' Main function that matches the original SecAct.activity.inference behavior,
#' including signature grouping and row expansion.
#'
#' @param inputProfile Gene expression matrix (genes x samples)
#' @param inputProfile_control Optional control expression matrix
#' @param is.differential Logical; if TRUE, inputProfile is already differential
#' @param is.paired Logical; if TRUE, perform paired differential calculation
#' @param is.singleSampleLevel Logical; if TRUE, calculate per-sample activity
#' @param sigMatrix Signature matrix path or "SecAct" for default
#' @param is.group.sig Logical; if TRUE, group similar signatures (default: TRUE)
#' @param is.group.cor Correlation threshold for grouping (default: 0.9)
#' @param lambda Ridge penalty factor (default: 5e+05)
#' @param nrand Number of permutations (default: 1000)
#' @param sigFilter Logical; if TRUE, filter signatures by available genes
#' @param ncores Number of CPU cores for parallel versions (default: auto)
#' @param method Which implementation to use: "legacy", "styp", "sttp", "mtyp", "mttp", "r".
#'   The "r" method uses pure R (no GSL required).
#' @param sort_genes Logical; if TRUE, sort common genes alphabetically before
#'   running ridge regression. This ensures reproducible results across different
#'   platforms but may differ from original gene order (default: FALSE).
#' @param verbose If TRUE, print progress messages
#'
#' @return Named list with beta, se, zscore, pvalue matrices (proteins x samples)
#'
#' @export
SecAct.activity.inference <- function(
    inputProfile,
    inputProfile_control = NULL,
    is.differential = FALSE,
    is.paired = FALSE,
    is.singleSampleLevel = FALSE,
    sigMatrix = "SecAct",
    is.group.sig = TRUE,
    is.group.cor = 0.9,
    lambda = 5e+05,
    nrand = 1000,
    sigFilter = FALSE,
    ncores = NULL,
    method = "sttp",
    sort_genes = FALSE,
    verbose = TRUE
) {
    if (verbose) {
        cat("SecAct Activity Inference\n")
        cat("==================================================\n")
        cat("  Input:", nrow(inputProfile), "genes ×", ncol(inputProfile), "samples\n")
    }
    
    # Compute differential expression if needed
    if (is.differential) {
        Y <- inputProfile
        if (ncol(Y) == 1) colnames(Y) <- "Change"
        if (verbose) cat("  Using pre-computed differential expression\n")
    } else {
        if (is.null(inputProfile_control)) {
            Y <- inputProfile - rowMeans(inputProfile)
            if (verbose) cat("  Differential: expr - rowMeans(expr)\n")
        } else {
            if (is.paired) {
                Y <- inputProfile - inputProfile_control[, colnames(inputProfile), drop = FALSE]
                if (verbose) cat("  Differential: paired comparison\n")
            } else {
                Y <- inputProfile - rowMeans(inputProfile_control)
                if (verbose) cat("  Differential: expr - mean(control)\n")
            }
            
            if (!is.singleSampleLevel) {
                Y <- matrix(rowMeans(Y), ncol = 1, 
                           dimnames = list(rownames(Y), "Change"))
                if (verbose) cat("  Aggregated to single column\n")
            }
        }
    }
    
    if (verbose) cat("  Expression matrix:", nrow(Y), "genes ×", ncol(Y), "samples\n")
    
    # Load signature matrix
    X <- .load_signature(sigMatrix)
    if (verbose) cat("  Loaded signature:", nrow(X), "genes ×", ncol(X), "proteins\n")
    
    # Filter signatures if requested
    if (sigFilter) {
        n_before <- ncol(X)
        X <- X[, colnames(X) %in% rownames(Y), drop = FALSE]
        if (verbose) cat("  sig_filter: kept", ncol(X), "/", n_before, "proteins\n")
    }
    
    # Group similar signatures if requested
    if (is.group.sig) {
        if (verbose) cat("  Grouping signatures (cor_threshold=", is.group.cor, ")...\n", sep = "")
        X <- .group_signatures(X, is.group.cor)
        if (verbose) cat("  Grouped into", ncol(X), "signature groups\n")
    }
    
    # Find overlapping genes
    olp <- intersect(rownames(Y), rownames(X))
    
    # Optionally sort genes alphabetically
    if (sort_genes) {
        olp <- sort(olp)
        if (verbose) cat("  Sorted genes alphabetically\n")
    }
    
    if (verbose) cat("  Common genes:", length(olp), "\n")
    
    if (length(olp) < 2) {
        stop("Too few overlapping genes between expression and signature matrices!")
    }
    
    # Prepare matrices
    X <- as.matrix(X[olp, , drop = FALSE])
    Y <- as.matrix(Y[olp, , drop = FALSE])
    
    X <- scale(X)
    Y <- scale(Y)
    
    X[is.na(X)] <- 0
    Y[is.na(Y)] <- 0
    
    if (verbose) {
        cat("  Running ridge regression (n_rand=", nrand, ")...\n", sep = "")
        cat("    X:", nrow(X), "genes ×", ncol(X), "features\n")
        cat("    Y:", nrow(Y), "genes ×", ncol(Y), "samples\n")
        cat("    lambda:", lambda, "\n")
        cat("    method:", method, "\n")
    }
    
    # Run inference using selected method
    ncores <- .default_ncores(ncores)
    
    start_time <- Sys.time()
    
    res <- switch(method,
        "legacy" = {
            n <- nrow(X)
            p <- ncol(X)
            m <- ncol(Y)
            len <- p * m
            
            # Legacy .C interface expects transposed matrices
            raw <- .C(
                "ridgeReg",
                as.double(t(X)),
                as.double(t(Y)),
                as.integer(n),
                as.integer(p),
                as.integer(m),
                as.double(lambda),
                as.double(nrand),
                beta   = double(len),
                se     = double(len),
                zscore = double(len),
                pvalue = double(len),
                PACKAGE = "RidgeR"
            )
            
            list(
                beta   = matrix(raw$beta, byrow = TRUE, ncol = m,
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, byrow = TRUE, ncol = m,
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, byrow = TRUE, ncol = m,
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, byrow = TRUE, ncol = m,
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        "styp" = {
            # .Call interface expects NON-transposed matrices
            raw <- .Call(
                "ridgeReg_styp_interface",
                X, Y,
                as.numeric(lambda),
                as.integer(nrand),
                PACKAGE = "RidgeR"
            )
            list(
                beta   = matrix(raw$beta, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        "sttp" = {
            # .Call interface expects NON-transposed matrices
            raw <- .Call(
                "ridgeReg_sttp_interface",
                X, Y,
                as.numeric(lambda),
                as.integer(nrand),
                PACKAGE = "RidgeR"
            )
            list(
                beta   = matrix(raw$beta, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        "mtyp" = {
            # .Call interface expects NON-transposed matrices
            raw <- .Call(
                "ridgeReg_mtyp_interface",
                X, Y,
                as.numeric(lambda),
                as.integer(nrand),
                as.integer(ncores),
                PACKAGE = "RidgeR"
            )
            list(
                beta   = matrix(raw$beta, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        "mttp" = {
            # .Call interface expects NON-transposed matrices
            raw <- .Call(
                "ridgeReg_mttp_interface",
                X, Y,
                as.numeric(lambda),
                as.integer(nrand),
                as.integer(ncores),
                PACKAGE = "RidgeR"
            )
            list(
                beta   = matrix(raw$beta, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, nrow = ncol(X), ncol = ncol(Y), byrow = TRUE,
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        "r" = {
            # Pure R implementation (no GSL required)
            n <- nrow(Y)
            p <- ncol(X)
            m <- ncol(Y)
            
            # Compute ridge regression via Cholesky decomposition
            A <- crossprod(X) + lambda * diag(p)
            R <- chol(A)
            beta <- backsolve(R, forwardsolve(t(R), crossprod(X, Y)))
            
            # Permutation testing
            aver <- NULL
            aver_sq <- NULL
            pvalue_mat <- NULL
            
            for (i in seq_len(nrand)) {
                set.seed(i)
                beta_rand <- backsolve(R, forwardsolve(t(R), crossprod(X, Y[sample.int(n), , drop = FALSE])))
                
                if (i == 1) {
                    aver <- beta_rand
                    aver_sq <- beta_rand^2
                    pvalue_mat <- (abs(beta_rand) >= abs(beta)) * 1.0
                } else {
                    aver <- aver + beta_rand
                    aver_sq <- aver_sq + beta_rand^2
                    pvalue_mat <- pvalue_mat + (abs(beta_rand) >= abs(beta))
                }
            }
            
            # Finalize statistics
            aver <- aver / nrand
            aver_sq <- aver_sq / nrand
            se <- sqrt(aver_sq - aver * aver)
            zscore <- (beta - aver) / se
            zscore[!is.finite(zscore)] <- 0
            pvalue_mat <- (pvalue_mat + 1) / (nrand + 1)
            
            # Set dimension names
            rownames(beta) <- colnames(X)
            colnames(beta) <- colnames(Y)
            rownames(se) <- colnames(X)
            colnames(se) <- colnames(Y)
            rownames(zscore) <- colnames(X)
            colnames(zscore) <- colnames(Y)
            rownames(pvalue_mat) <- colnames(X)
            colnames(pvalue_mat) <- colnames(Y)
            
            list(beta = beta, se = se, zscore = zscore, pvalue = pvalue_mat)
        },
        stop("Unknown method: ", method, ". Valid methods: legacy, styp, sttp, mtyp, mttp, r")
    )
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (verbose) cat("  Ridge regression completed in", round(elapsed, 2), "seconds\n")
    
    # Expand grouped signatures back to individual rows
    if (is.group.sig) {
        if (verbose) cat("  Expanding grouped signatures...\n")
        res$beta   <- .expand_rows(res$beta)
        res$se     <- .expand_rows(res$se)
        res$zscore <- .expand_rows(res$zscore)
        res$pvalue <- .expand_rows(res$pvalue)
        
        # Sort by row name
        row_order <- sort(rownames(res$beta))
        res$beta   <- res$beta[row_order, , drop = FALSE]
        res$se     <- res$se[row_order, , drop = FALSE]
        res$zscore <- res$zscore[row_order, , drop = FALSE]
        res$pvalue <- res$pvalue[row_order, , drop = FALSE]
    }
    
    if (verbose) {
        cat("  Result shape:", nrow(res$beta), "proteins ×", ncol(res$beta), "samples\n")
        cat("==================================================\n")
    }
    
    res
}


# ==============================================================================
# SECTION 7: TESTING & COMPARISON UTILITIES
# ==============================================================================

#' Compare All Implementation Methods
#'
#' Utility function to verify that all implementations produce identical
#' (within tolerance) results and to compare their execution times.
#'
#' @inheritParams SecAct.inference.gsl.mtyp
#' @param tol Tolerance for numerical comparison (default: 1e-10)
#'
#' @return Invisibly returns a named list containing results from all methods
#'
#' @export
SecAct.compare.methods <- function(Y,
                                    SigMat = "SecAct",
                                    lambda = 5e+05,
                                    nrand = 100,
                                    ncores = NULL,
                                    tol = 1e-10) {

    # Helper to time expressions
    time_it <- function(expr) {
        system.time(expr)["elapsed"]
    }

    # Run all implementations
    cat("Running gsl.legacy (.C legacy interface)...\n")
    t_legacy <- time_it({
        res_legacy <- SecAct.inference.gsl.legacy(Y, SigMat, lambda, nrand)
    })

    cat("Running gsl.styp (single-threaded, Y-permutation)...\n")
    t_styp <- time_it({
        res_styp <- SecAct.inference.gsl.styp(Y, SigMat, lambda, nrand)
    })

    cat("Running gsl.sttp (single-threaded, T-permutation)...\n")
    t_sttp <- time_it({
        res_sttp <- SecAct.inference.gsl.sttp(Y, SigMat, lambda, nrand)
    })

    cat("Running gsl.mtyp (multi-threaded, Y-permutation)...\n")
    t_mtyp <- time_it({
        res_mtyp <- SecAct.inference.gsl.mtyp(Y, SigMat, lambda, nrand, ncores)
    })

    cat("Running gsl.mttp (multi-threaded, T-permutation)...\n")
    t_mttp <- time_it({
        res_mttp <- SecAct.inference.gsl.mttp(Y, SigMat, lambda, nrand, ncores)
    })

    cat("Running pure R implementation...\n")
    t_r <- time_it({
        res_r <- SecAct.inference.r(Y, SigMat, lambda, nrand)
    })

    # Helper to check equality within tolerance
    check_equal <- function(a, b) {
        max_diff <- max(abs(a - b), na.rm = TRUE)
        max_diff < tol
    }

    # Report consistency
    cat("\n")
    cat("=== Consistency Check (z-score comparison) ===\n")
    cat(sprintf("  legacy vs styp : %s\n",
                ifelse(check_equal(res_legacy$zscore, res_styp$zscore), "PASS", "FAIL")))
    cat(sprintf("  styp   vs mtyp : %s\n",
                ifelse(check_equal(res_styp$zscore, res_mtyp$zscore), "PASS", "FAIL")))
    cat(sprintf("  sttp   vs mttp : %s\n",
                ifelse(check_equal(res_sttp$zscore, res_mttp$zscore), "PASS", "FAIL")))
    cat(sprintf("  mtyp   vs mttp : %s\n",
                ifelse(check_equal(res_mtyp$zscore, res_mttp$zscore), "PASS", "FAIL")))

    # Report timing
    cat("\n")
    cat("=== Timing Summary ===\n")
    cat(sprintf("  gsl.legacy : %6.2f seconds\n", t_legacy))
    cat(sprintf("  gsl.styp   : %6.2f seconds\n", t_styp))
    cat(sprintf("  gsl.sttp   : %6.2f seconds\n", t_sttp))
    cat(sprintf("  gsl.mtyp   : %6.2f seconds\n", t_mtyp))
    cat(sprintf("  gsl.mttp   : %6.2f seconds\n", t_mttp))
    cat(sprintf("  pure R     : %6.2f seconds\n", t_r))

    # Return all results invisibly
    invisible(list(
        legacy = res_legacy,
        styp   = res_styp,
        sttp   = res_sttp,
        mtyp   = res_mtyp,
        mttp   = res_mttp,
        r      = res_r
    ))
}

#' Compare RidgeR with SecAct Output
#'
#' Utility function to compare SecAct.activity.inference results across
#' different methods and optionally with external SecAct package results.
#'
#' @param inputProfile Gene expression matrix (genes x samples)
#' @param is.differential Logical; if TRUE, inputProfile is already differential
#' @param sigMatrix Signature matrix path or "SecAct" for default
#' @param is.group.sig Logical; if TRUE, group similar signatures
#' @param is.group.cor Correlation threshold for grouping
#' @param lambda Ridge penalty factor
#' @param nrand Number of permutations
#' @param ncores Number of CPU cores
#'
#' @return Invisibly returns list of results from all methods
#'
#' @examples
#' \dontrun{
#' dataPath <- file.path(system.file(package = "RidgeR"), "extdata")
#' expr.diff <- read.table(paste0(dataPath, "/Ly86-Fc_vs_Vehicle_logFC.txt"))
#' results <- SecAct.compare.activity(expr.diff, is.differential = TRUE)
#' }
#'
#' @export
SecAct.compare.activity <- function(
    inputProfile,
    is.differential = FALSE,
    sigMatrix = "SecAct",
    is.group.sig = TRUE,
    is.group.cor = 0.9,
    lambda = 5e+05,
    nrand = 1000,
    ncores = NULL
) {
    # Helper to time expressions
    time_it <- function(expr) {
        system.time(expr)["elapsed"]
    }
    
    cat("Running SecAct.activity.inference with different methods...\n\n")
    
    # Run all implementations
    cat("Method: legacy (.C interface)...\n")
    t_legacy <- time_it({
        res_legacy <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, method = "legacy"
        )
    })
    
    cat("Method: styp (single-threaded, Y-permutation)...\n")
    t_styp <- time_it({
        res_styp <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, method = "styp"
        )
    })
    
    cat("Method: sttp (single-threaded, T-permutation)...\n")
    t_sttp <- time_it({
        res_sttp <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, method = "sttp"
        )
    })
    
    cat("Method: mtyp (multi-threaded, Y-permutation)...\n")
    t_mtyp <- time_it({
        res_mtyp <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, ncores = ncores, method = "mtyp"
        )
    })
    
    cat("Method: mttp (multi-threaded, T-permutation)...\n")
    t_mttp <- time_it({
        res_mttp <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, ncores = ncores, method = "mttp"
        )
    })
    
    cat("Method: r (pure R, no GSL)...\n")
    t_r <- time_it({
        res_r <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, method = "r"
        )
    })
    
    # Report results
    cat("\n=== Z-score Comparison (first 6 rows) ===\n")
    cat("\nLegacy:\n")
    print(head(res_legacy$zscore))
    cat("\nMTYP:\n")
    print(head(res_mtyp$zscore))
    cat("\nPure R:\n")
    print(head(res_r$zscore))
    
    # Check consistency
    tol <- 1e-10
    check_equal <- function(a, b) {
        max(abs(a - b), na.rm = TRUE) < tol
    }
    
    cat("\n=== Consistency Check ===\n")
    cat(sprintf("  legacy vs styp : %s\n",
                ifelse(check_equal(res_legacy$zscore, res_styp$zscore), "PASS", "FAIL")))
    cat(sprintf("  styp   vs mtyp : %s\n",
                ifelse(check_equal(res_styp$zscore, res_mtyp$zscore), "PASS", "FAIL")))
    cat(sprintf("  sttp   vs mttp : %s\n",
                ifelse(check_equal(res_sttp$zscore, res_mttp$zscore), "PASS", "FAIL")))
    cat(sprintf("  mtyp   vs mttp : %s\n",
                ifelse(check_equal(res_mtyp$zscore, res_mttp$zscore), "PASS", "FAIL")))
    
    cat("\n=== Timing Summary ===\n")
    cat(sprintf("  legacy : %6.2f seconds\n", t_legacy))
    cat(sprintf("  styp   : %6.2f seconds\n", t_styp))
    cat(sprintf("  sttp   : %6.2f seconds\n", t_sttp))
    cat(sprintf("  mtyp   : %6.2f seconds\n", t_mtyp))
    cat(sprintf("  mttp   : %6.2f seconds\n", t_mttp))
    cat(sprintf("  r      : %6.2f seconds\n", t_r))
    
    invisible(list(
        legacy = res_legacy,
        styp   = res_styp,
        sttp   = res_sttp,
        mtyp   = res_mtyp,
        mttp   = res_mttp,
        r      = res_r
    ))
}

# ==============================================================================
# SECTION 8: ST DATA HELPER FUNCTIONS
# ==============================================================================

#' Transfer gene symbols to standard format
#'
#' Converts gene symbols to uppercase and handles common naming variations.
#'
#' @param symbols Character vector of gene symbols
#' @return Character vector of standardized gene symbols
#' @keywords internal
.transfer_symbol <- function(symbols) {
    # Convert to uppercase for consistency
    symbols <- toupper(symbols)
    
    # Remove version numbers (e.g., "GENE.1" -> "GENE")
    symbols <- gsub("\\.\\d+$", "", symbols)
    
    symbols
}


#' Remove duplicate rows by keeping the one with highest mean
#'
#' When multiple rows have the same gene name, keeps the row with
#' the highest mean expression value.
#'
#' @param mat Expression matrix with gene names as rownames
#' @return Matrix with unique row names
#' @keywords internal
.rm_duplicates <- function(mat) {
    if (!any(duplicated(rownames(mat)))) {
        return(mat)
    }
    
    # Calculate row means
    if (inherits(mat, "dgCMatrix") || inherits(mat, "sparseMatrix")) {
        row_means <- Matrix::rowMeans(mat)
    } else {
        row_means <- rowMeans(mat)
    }
    
    # For each unique gene, keep the row with highest mean
    unique_genes <- unique(rownames(mat))
    keep_idx <- sapply(unique_genes, function(gene) {
        idx <- which(rownames(mat) == gene)
        if (length(idx) == 1) return(idx)
        idx[which.max(row_means[idx])]
    })
    
    mat[keep_idx, , drop = FALSE]
}


#' Sweep operation for sparse matrices
#'
#' Applies a sweep operation to sparse matrices without converting to dense.
#'
#' @param x Sparse matrix
#' @param MARGIN 1 for rows, 2 for columns
#' @param STATS Vector of statistics to apply
#' @param FUN Function to apply (typically "/" or "-")
#' @return Sparse matrix with operation applied
#' @keywords internal
.sweep_sparse <- function(x, MARGIN, STATS, FUN = "/") {
    if (!inherits(x, "dgCMatrix") && !inherits(x, "sparseMatrix")) {
        return(sweep(x, MARGIN, STATS, FUN))
    }
    
    if (MARGIN == 2) {
        if (FUN == "/") {
            # Column-wise division
            x@x <- x@x / rep(STATS, diff(x@p))
        } else if (FUN == "-") {
            # Column-wise subtraction - need to handle differently
            x <- as(x, "dgCMatrix")
            for (j in seq_along(STATS)) {
                idx <- (x@p[j] + 1):x@p[j + 1]
                if (length(idx) > 0) {
                    x@x[idx] <- x@x[idx] - STATS[j]
                }
            }
        }
    } else if (MARGIN == 1) {
        if (FUN == "/") {
            # Row-wise division - less efficient for CSC format
            x <- t(.sweep_sparse(t(x), 2, STATS, FUN))
        } else if (FUN == "-") {
            x <- t(.sweep_sparse(t(x), 2, STATS, FUN))
        }
    }
    
    x
}


# ==============================================================================
# SECTION 9: SPATIAL TRANSCRIPTOMICS ACTIVITY INFERENCE
# ==============================================================================

#' Spot Activity Inference from Spatial Transcriptomics Data
#'
#' Calculate secreted protein signaling activity of spots from spatial
#' transcriptomics data. Supports SpaCET objects and raw count matrices.
#'
#' @param inputProfile Either a SpaCET object or a gene expression count matrix
#'   (genes x spots). If a matrix, rownames should be gene symbols.
#' @param inputProfile_control Optional control expression data.
#' @param scale.factor Sets the scale factor for spot-level normalization
#'   (default: 1e+05, i.e., counts per 100K).
#' @param sigMatrix Secreted protein signature matrix path or "SecAct" for default.
#' @param is.group.sig Logical indicating whether to group similar signatures
#'   (default: TRUE).
#' @param is.group.cor Correlation cutoff for signature grouping (default: 0.9).
#' @param lambda Ridge penalty factor (default: 5e+05).
#' @param nrand Number of permutations (default: 1000).
#' @param sigFilter Logical indicating whether to filter signatures by available
#'   genes (default: FALSE).
#' @param ncores Number of CPU cores for parallel processing (default: auto).
#' @param method Which implementation to use (default: "sttp").
#' @param sort_genes Logical; if TRUE, sort common genes alphabetically (default: FALSE).
#' @param return.SpaCET If inputProfile is a SpaCET object, whether to return
#'   the modified SpaCET object (TRUE) or just the results list (FALSE).
#'
#' @return SpaCET object with results or named list with beta, se, zscore, pvalue matrices.
#'
#' @export
SecAct.activity.inference.ST <- function(
    inputProfile,
    inputProfile_control = NULL,
    scale.factor = 1e+05,
    sigMatrix = "SecAct",
    is.group.sig = TRUE,
    is.group.cor = 0.9,
    lambda = 5e+05,
    nrand = 1000,
    sigFilter = FALSE,
    ncores = NULL,
    method = "sttp",
    sort_genes = FALSE,
    return.SpaCET = TRUE
) {
    # Check if input is a SpaCET object
    is_spacet <- inherits(inputProfile, "SpaCET")
    
    if (is_spacet) {
        # Extract count matrix from SpaCET object
        expr <- inputProfile@input$counts
        expr <- expr[Matrix::rowSums(expr) > 0, ]
        rownames(expr) <- .transfer_symbol(rownames(expr))
        expr <- .rm_duplicates(expr)
    } else {
        # Assume inputProfile is a count matrix
        if (inherits(inputProfile, "sparseMatrix") || inherits(inputProfile, "dgCMatrix")) {
            expr <- inputProfile
            expr <- expr[Matrix::rowSums(expr) > 0, ]
        } else {
            expr <- as.matrix(inputProfile)
            expr <- expr[rowSums(expr) > 0, ]
        }
        rownames(expr) <- .transfer_symbol(rownames(expr))
        expr <- .rm_duplicates(expr)
    }
    
    # Normalize to counts per scale.factor
    if (inherits(expr, "sparseMatrix") || inherits(expr, "dgCMatrix")) {
        stats <- Matrix::colSums(expr)
        expr <- .sweep_sparse(expr, 2, stats, "/")
        expr@x <- expr@x * scale.factor
        # Log2 transform
        expr@x <- log2(expr@x + 1)
    } else {
        stats <- colSums(expr)
        expr <- sweep(expr, 2, stats, "/") * scale.factor
        # Log2 transform
        expr <- log2(expr + 1)
    }
    
    # Compute differential expression
    if (is.null(inputProfile_control)) {
        # Use mean of inputProfile as control
        if (inherits(expr, "sparseMatrix") || inherits(expr, "dgCMatrix")) {
            expr.diff <- expr - Matrix::rowMeans(expr)
        } else {
            expr.diff <- expr - rowMeans(expr)
        }
    } else {
        # Process control expression
        if (inherits(inputProfile_control, "SpaCET")) {
            expr_control <- inputProfile_control@input$counts
            expr_control <- expr_control[Matrix::rowSums(expr_control) > 0, ]
            rownames(expr_control) <- .transfer_symbol(rownames(expr_control))
            expr_control <- .rm_duplicates(expr_control)
        } else {
            if (inherits(inputProfile_control, "sparseMatrix") || 
                inherits(inputProfile_control, "dgCMatrix")) {
                expr_control <- inputProfile_control
                expr_control <- expr_control[Matrix::rowSums(expr_control) > 0, ]
            } else {
                expr_control <- as.matrix(inputProfile_control)
                expr_control <- expr_control[rowSums(expr_control) > 0, ]
            }
            rownames(expr_control) <- .transfer_symbol(rownames(expr_control))
            expr_control <- .rm_duplicates(expr_control)
        }
        
        # Normalize control
        if (inherits(expr_control, "sparseMatrix") || inherits(expr_control, "dgCMatrix")) {
            stats <- Matrix::colSums(expr_control)
            expr_control <- .sweep_sparse(expr_control, 2, stats, "/")
            expr_control@x <- expr_control@x * scale.factor
            expr_control@x <- log2(expr_control@x + 1)
        } else {
            stats <- colSums(expr_control)
            expr_control <- sweep(expr_control, 2, stats, "/") * scale.factor
            expr_control <- log2(expr_control + 1)
        }
        
        # Find overlapping genes and compute difference
        olp <- intersect(rownames(expr), rownames(expr_control))
        if (inherits(expr, "sparseMatrix") || inherits(expr, "dgCMatrix")) {
            expr.diff <- expr[olp, ] - Matrix::rowMeans(expr_control[olp, ])
        } else {
            expr.diff <- expr[olp, ] - rowMeans(expr_control[olp, ])
        }
    }
    
    # Convert sparse to dense for ridge regression (if small enough)
    if (inherits(expr.diff, "sparseMatrix") || inherits(expr.diff, "dgCMatrix")) {
        expr.diff <- as.matrix(expr.diff)
    }
    
    # Run activity inference
    res <- SecAct.activity.inference(
        inputProfile = expr.diff,
        is.differential = TRUE,
        sigMatrix = sigMatrix,
        is.group.sig = is.group.sig,
        is.group.cor = is.group.cor,
        lambda = lambda,
        nrand = nrand,
        sigFilter = sigFilter,
        ncores = ncores,
        method = method,
        sort_genes = sort_genes
    )
    
    # Return results
    if (is_spacet && return.SpaCET) {
        # Store results in SpaCET object
        if (is.null(inputProfile@results)) {
            inputProfile@results <- list()
        }
        if (is.null(inputProfile@results$SecAct_output)) {
            inputProfile@results$SecAct_output <- list()
        }
        inputProfile@results$SecAct_output$SecretedProteinActivity <- res
        return(inputProfile)
    } else {
        return(res)
    }
}


#' Cell State Activity Inference from Single Cell RNA-seq Data
#'
#' Calculate secreted protein signaling activity of cell states from
#' single cell RNA-Sequencing data. Supports Seurat objects.
#'
#' @param inputProfile A Seurat object.
#' @param cellType_meta Column name in meta.data containing cell-type annotations.
#' @param sigMatrix Secreted protein signature matrix path or "SecAct" for default.
#' @param is.singleCellLevel Logical indicating whether to calculate activity
#'   for each single cell (TRUE) or aggregate by cell type (FALSE, default).
#' @param is.group.sig Logical indicating whether to group similar signatures.
#' @param is.group.cor Correlation cutoff for signature grouping.
#' @param lambda Ridge penalty factor.
#' @param nrand Number of permutations.
#' @param sigFilter Logical indicating whether to filter signatures.
#' @param ncores Number of CPU cores for parallel processing.
#' @param method Which implementation to use (default: "sttp").
#' @param sort_genes Logical; if TRUE, sort common genes alphabetically (default: FALSE).
#' @param return.Seurat Whether to return the modified Seurat object (TRUE)
#'   or just the results list (FALSE).
#'
#' @return Seurat object with results or named list with beta, se, zscore, pvalue matrices.
#'
#' @export
SecAct.activity.inference.scRNAseq <- function(
    inputProfile,
    cellType_meta,
    sigMatrix = "SecAct",
    is.singleCellLevel = FALSE,
    is.group.sig = TRUE,
    is.group.cor = 0.9,
    lambda = 5e+05,
    nrand = 1000,
    sigFilter = FALSE,
    ncores = NULL,
    method = "sttp",
    sort_genes = FALSE,
    return.Seurat = TRUE
) {
    # Check if input is a Seurat object
    if (!inherits(inputProfile, "Seurat")) {
        stop("Please input a Seurat object.")
    }
    
    # Extract count matrix - handle Seurat v5 (Assay5) vs older versions
    if (inherits(inputProfile@assays$RNA, "Assay5")) {
        counts <- inputProfile@assays$RNA@layers$counts
        colnames(counts) <- rownames(inputProfile@assays$RNA@cells)
        rownames(counts) <- rownames(inputProfile@assays$RNA@features)
    } else {
        counts <- inputProfile@assays$RNA@counts
    }
    
    # Standardize gene symbols
    rownames(counts) <- .transfer_symbol(rownames(counts))
    counts <- .rm_duplicates(counts)
    
    if (!is.singleCellLevel) {
        # Generate pseudo-bulk by cell type
        cellType_vec <- inputProfile@meta.data[, cellType_meta]
        
        expr <- data.frame(row.names = rownames(counts))
        for (cellType in sort(unique(cellType_vec))) {
            cells_idx <- which(cellType_vec == cellType)
            if (inherits(counts, "sparseMatrix") || inherits(counts, "dgCMatrix")) {
                expr[, cellType] <- Matrix::rowSums(counts[, cells_idx, drop = FALSE])
            } else {
                expr[, cellType] <- rowSums(counts[, cells_idx, drop = FALSE])
            }
        }
        
        # Normalize to TPM (counts per million)
        if (inherits(expr, "sparseMatrix") || inherits(expr, "dgCMatrix")) {
            expr <- sweep(as.matrix(expr), 2, colSums(expr), "/") * 1e6
        } else {
            expr <- sweep(expr, 2, colSums(expr), "/") * 1e6
        }
    } else {
        # Single cell level - normalize per cell
        if (inherits(counts, "sparseMatrix") || inherits(counts, "dgCMatrix")) {
            stats <- Matrix::colSums(counts)
            expr <- .sweep_sparse(counts, 2, stats, "/")
            expr@x <- expr@x * 1e5
            expr <- as.matrix(expr)
        } else {
            expr <- sweep(counts, 2, colSums(counts), "/") * 1e5
        }
    }
    rm(counts)
    gc()
    
    # Log2 transform
    expr <- log2(expr + 1)
    
    # Compute differential expression (vs mean)
    expr.diff <- expr - rowMeans(expr)
    rm(expr)
    gc()
    
    # Run activity inference
    res <- SecAct.activity.inference(
        inputProfile = expr.diff,
        is.differential = TRUE,
        sigMatrix = sigMatrix,
        is.group.sig = is.group.sig,
        is.group.cor = is.group.cor,
        lambda = lambda,
        nrand = nrand,
        sigFilter = sigFilter,
        ncores = ncores,
        method = method,
        sort_genes = sort_genes
    )
    
    # Return results
    if (return.Seurat) {
        if (is.null(inputProfile@misc)) {
            inputProfile@misc <- list()
        }
        if (is.null(inputProfile@misc$SecAct_output)) {
            inputProfile@misc$SecAct_output <- list()
        }
        inputProfile@misc$SecAct_output$SecretedProteinActivity <- res
        return(inputProfile)
    } else {
        return(res)
    }
}


# ==============================================================================
# SECTION 10: H5AD I/O UTILITIES
# ==============================================================================

#' Write SecAct Results to H5AD Format (Python-compatible)
#'
#' Exports SecAct results (beta, se, zscore, pvalue) from SpaCET, Seurat objects,
#' or a list to H5AD format for interoperability with Python/scanpy.
#' 
#' Uses reticulate + anndata for fully valid H5AD output.
#'
#' @param obj A SpaCET object, Seurat object, or a list containing SecAct results.
#'   If a list, it should contain: beta, se, zscore, pvalue matrices (proteins x cells).
#' @param output_file Character string. Output file path (default: "SecAct_results.h5ad").
#' @param compression Character string. Compression method: "gzip" or NULL (default: "gzip").
#'
#' @return Invisibly returns the output file path.
#'
#' @details
#' The function creates an anndata-compatible H5AD file with:
#' \itemize{
#'   \item \code{X}: Beta coefficients matrix (cells x proteins)
#'   \item \code{obsm/se}: Standard errors (cells x proteins)
#'   \item \code{obsm/zscore}: Z-scores (cells x proteins)
#'   \item \code{obsm/pvalue}: P-values (cells x proteins)
#'   \item \code{uns/source}: Source information
#' }
#'
#' @section Python Usage:
#' \preformatted{
#' import anndata
#' adata = anndata.read_h5ad("SecAct_results.h5ad")
#' beta   = adata.X                    # (n_cells, n_proteins)
#' se     = adata.obsm['se']           # (n_cells, n_proteins)
#' zscore = adata.obsm['zscore']       # (n_cells, n_proteins)
#' pvalue = adata.obsm['pvalue']       # (n_cells, n_proteins)
#' protein_names = adata.var_names     # protein names
#' cell_names    = adata.obs_names     # cell/sample names
#' }
#'
#' @section Requirements:
#' \preformatted{
#' # Install reticulate
#' install.packages("reticulate")
#' 
#' # Install anndata in Python
#' reticulate::py_install("anndata", pip = TRUE)
#' }
#'
#' @examples
#' \dontrun{
#' # From SpaCET object
#' write_secact_to_h5ad(SpaCET_obj, "spacet_results.h5ad")
#'
#' # From Seurat object
#' write_secact_to_h5ad(seurat_obj, "seurat_results.h5ad")
#'
#' # From list (direct SecAct output)
#' res <- list(beta = beta_mat, se = se_mat, zscore = zscore_mat, pvalue = pvalue_mat)
#' write_secact_to_h5ad(res, "results.h5ad")
#' }
#'
#' @export
write_secact_to_h5ad <- function(obj, output_file = "SecAct_results.h5ad", compression = "gzip") {

    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required. Install with: install.packages('reticulate')")
    }

    cat("Saving SecAct results to H5AD (via Python anndata)\n")
    cat("Output:", output_file, "\n\n")

    ## ----------------------------------------------------------------
    ## 1. Detect object type & extract SecAct results
    ## ----------------------------------------------------------------
    if (inherits(obj, "SpaCET")) {
        cat("Detected SpaCET object\n")
        if (is.null(obj@results$SecAct_output$SecretedProteinActivity)) {
            stop("No SecAct results found in SpaCET object. Run SecAct inference first.")
        }
        res <- obj@results$SecAct_output$SecretedProteinActivity
        source <- "SpaCET / SecAct"

    } else if (inherits(obj, "Seurat")) {
        cat("Detected Seurat object\n")
        if (is.null(obj@misc$SecAct_output$SecretedProteinActivity)) {
            stop("No SecAct results found in Seurat object. Run SecAct inference first.")
        }
        res <- obj@misc$SecAct_output$SecretedProteinActivity
        source <- "Seurat / SecAct"

    } else if (is.list(obj)) {
        cat("Detected list object\n")
        # Validate list has required components
        required <- c("beta", "se", "zscore", "pvalue")
        missing <- setdiff(required, names(obj))
        if (length(missing) > 0) {
            stop("List is missing required components: ", paste(missing, collapse = ", "),
                 "\nRequired: beta, se, zscore, pvalue")
        }
        res <- obj
        source <- "List / SecAct"

    } else {
        stop("Unsupported object type. Must be SpaCET, Seurat, or a list with beta/se/zscore/pvalue.")
    }

    ## ----------------------------------------------------------------
    ## 2. Extract matrices (R format: proteins x cells)
    ##    Handle vector, 1-column matrix, and multi-column matrix inputs
    ## ----------------------------------------------------------------
    beta_r   <- res$beta
    se_r     <- res$se
    zscore_r <- res$zscore
    pvalue_r <- res$pvalue

    # Handle vector input (single sample case - no dim attribute)
    if (is.vector(beta_r) && is.null(dim(beta_r))) {
        cat("Detected single-sample data (vectors). Converting to matrices...\n")
        
        # Get protein names from vector names
        protein_names <- names(beta_r)
        if (is.null(protein_names)) {
            protein_names <- paste0("protein_", seq_along(beta_r))
            warning("No protein names found. Using auto-generated names.")
        }
        
        # Convert to column matrix (proteins x 1 sample)
        beta_r   <- matrix(beta_r, ncol = 1, dimnames = list(protein_names, "sample_1"))
        se_r     <- matrix(se_r, ncol = 1, dimnames = list(protein_names, "sample_1"))
        zscore_r <- matrix(zscore_r, ncol = 1, dimnames = list(protein_names, "sample_1"))
        pvalue_r <- matrix(pvalue_r, ncol = 1, dimnames = list(protein_names, "sample_1"))
    } else {
        # Convert to matrix
        beta_r   <- as.matrix(beta_r)
        se_r     <- as.matrix(se_r)
        zscore_r <- as.matrix(zscore_r)
        pvalue_r <- as.matrix(pvalue_r)
        
        # Handle 1-column matrix with no colnames (single sample stored as matrix)
        if (ncol(beta_r) == 1 && is.null(colnames(beta_r))) {
            cat("Detected single-sample data (1-column matrix). Setting sample name...\n")
            colnames(beta_r)   <- "sample_1"
            colnames(se_r)     <- "sample_1"
            colnames(zscore_r) <- "sample_1"
            colnames(pvalue_r) <- "sample_1"
        }
    }

    stopifnot(
        "Dimension mismatch between beta and se" = identical(dim(beta_r), dim(se_r)),
        "Dimension mismatch between beta and zscore" = identical(dim(beta_r), dim(zscore_r)),
        "Dimension mismatch between beta and pvalue" = identical(dim(beta_r), dim(pvalue_r))
    )

    # Get names from matrix dimensions
    protein_names <- rownames(beta_r)
    cell_names <- colnames(beta_r)
    
    # Generate names if missing
    if (is.null(protein_names)) {
        protein_names <- paste0("protein_", seq_len(nrow(beta_r)))
        warning("No protein names found in rownames(beta). Using auto-generated names.")
    }
    if (is.null(cell_names)) {
        cell_names <- paste0("sample_", seq_len(ncol(beta_r)))
        warning("No sample/cell names found in colnames(beta). Using auto-generated names.")
    }
    
    n_proteins <- nrow(beta_r)
    n_cells <- ncol(beta_r)

    cat("Cells   :", n_cells, "\n")
    cat("Proteins:", n_proteins, "\n\n")

    ## ----------------------------------------------------------------
    ## 3. Transpose to Python format (cells x proteins)
    ## ----------------------------------------------------------------
    cat("Transposing to Python format (cells x proteins)...\n")
    
    beta   <- t(beta_r)      # cells x proteins
    se     <- t(se_r)        # cells x proteins
    zscore <- t(zscore_r)    # cells x proteins
    pvalue <- t(pvalue_r)    # cells x proteins

    cat("Output shape: (", n_cells, " x ", n_proteins, ") = (cells x proteins)\n\n", sep = "")

    ## ----------------------------------------------------------------
    ## 4. Write H5AD using anndata (CORRECT way)
    ## ----------------------------------------------------------------
    cat("Writing H5AD via Python anndata...\n")
    
    # Remove existing file if present
    if (file.exists(output_file)) file.remove(output_file)
    
    # Use Python directly to avoid reticulate type conversion issues
    # Import modules
    reticulate::py_run_string("import anndata")
    reticulate::py_run_string("import numpy as np")
    
    # Get the Python main module to pass data
    py <- reticulate::import("__main__")
    
    # Pass matrices to Python
    py$X_data <- beta
    py$se_data <- se
    py$zscore_data <- zscore
    py$pvalue_data <- pvalue
    py$cell_names <- as.list(cell_names)
    py$protein_names <- as.list(protein_names)
    py$output_file <- output_file
    py$source_info <- source
    py$use_compression <- !is.null(compression) && compression == "gzip"
    
    # Run Python code to create and save AnnData
    reticulate::py_run_string("
# Convert to numpy arrays and ensure 2D
X_np = np.atleast_2d(np.array(X_data, dtype=np.float64))
se_np = np.atleast_2d(np.array(se_data, dtype=np.float64))
zscore_np = np.atleast_2d(np.array(zscore_data, dtype=np.float64))
pvalue_np = np.atleast_2d(np.array(pvalue_data, dtype=np.float64))

print(f'  X shape: {X_np.shape}')
print(f'  se shape: {se_np.shape}')
print(f'  zscore shape: {zscore_np.shape}')
print(f'  pvalue shape: {pvalue_np.shape}')

# Create AnnData
adata = anndata.AnnData(X=X_np)

# Set names
adata.obs_names = cell_names
adata.var_names = protein_names

# Add obsm
adata.obsm['se'] = se_np
adata.obsm['zscore'] = zscore_np
adata.obsm['pvalue'] = pvalue_np

print(f'  obsm keys: {list(adata.obsm.keys())}')

# Add metadata
adata.uns['source'] = source_info
adata.uns['data_type'] = 'SecAct results'
adata.uns['description'] = 'X=beta, obsm={se,zscore,pvalue}'

# Write file
if use_compression:
    adata.write_h5ad(output_file, compression='gzip')
else:
    adata.write_h5ad(output_file)
")

    cat("\nDONE\n")
    cat("File size:", sprintf("%.2f MB", file.info(output_file)$size / 1e6), "\n")

    cat("\nPython usage:\n")
    cat("  import anndata\n")
    cat("  adata = anndata.read_h5ad('", basename(output_file), "')\n", sep = "")
    cat("  beta   = adata.X                    # shape: (", n_cells, ", ", n_proteins, ")\n", sep = "")
    cat("  se     = adata.obsm['se']           # shape: (", n_cells, ", ", n_proteins, ")\n", sep = "")
    cat("  zscore = adata.obsm['zscore']       # shape: (", n_cells, ", ", n_proteins, ")\n", sep = "")
    cat("  pvalue = adata.obsm['pvalue']       # shape: (", n_cells, ", ", n_proteins, ")\n", sep = "")
    cat("  proteins = list(adata.var_names)   # ", n_proteins, " protein names\n", sep = "")
    cat("  cells    = list(adata.obs_names)   # ", n_cells, " cell names\n", sep = "")

    invisible(output_file)
}


#' Read H5AD File into Seurat or SpaCET Object
#'
#' Read an AnnData (.h5ad) file and convert it to a Seurat object
#' (no spatial data) or SpaCET object (if spatial data available).
#'
#' Uses reticulate + anndata for proper H5AD reading.
#'
#' @param h5ad_file Path to .h5ad file
#'
#' @return A Seurat or SpaCET object with SecAct results in misc/results slot.
#'
#' @details
#' This function reads H5AD files created by \code{write_secact_to_h5ad()} or
#' Python SecActPy and reconstructs the SecAct result matrices.
#'
#' @examples
#' \dontrun{
#' # Read H5AD file
#' obj <- read_h5ad_to_secact("SecAct_results.h5ad")
#'
#' # Access results (always returns R format: proteins x cells)
#' res <- obj@misc$SecAct_output$SecretedProteinActivity
#' zscore <- res$zscore
#' }
#'
#' @export
read_h5ad_to_secact <- function(h5ad_file) {

    if (!file.exists(h5ad_file)) {
        stop("h5ad_file does not exist: ", h5ad_file)
    }

    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required. Install with: install.packages('reticulate')")
    }
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required. Install with: install.packages('Seurat')")
    }

    cat("Reading H5AD file:", h5ad_file, "\n")

    ## ----------------------------------------------------------------
    ## 1. Read H5AD using anndata
    ## ----------------------------------------------------------------
    anndata <- reticulate::import("anndata", delay_load = FALSE)
    np <- reticulate::import("numpy", delay_load = FALSE)
    
    adata <- anndata$read_h5ad(h5ad_file)
    
    cat("Shape:", adata$shape[[1]], "x", adata$shape[[2]], "\n")
    cat("obs (rows):", length(adata$obs_names), "\n")
    cat("var (cols):", length(adata$var_names), "\n")

    ## ----------------------------------------------------------------
    ## 2. Extract matrices (Python format: cells x proteins)
    ## ----------------------------------------------------------------
    # Get X as dense numpy array, then convert to R
    X_py <- adata$X
    
    # Handle sparse matrices
    scipy_sparse <- reticulate::import("scipy.sparse", delay_load = FALSE)
    if (scipy_sparse$issparse(X_py)) {
        X_py <- X_py$toarray()
    }
    beta_py <- reticulate::py_to_r(np$array(X_py))
    
    # Get obs/var names
    cell_names <- reticulate::py_to_r(adata$obs_names$tolist())
    protein_names <- reticulate::py_to_r(adata$var_names$tolist())
    
    # Get obsm matrices
    get_obsm_matrix <- function(name) {
        if (name %in% names(adata$obsm)) {
            mat <- adata$obsm[[name]]
            if (scipy_sparse$issparse(mat)) {
                mat <- mat$toarray()
            }
            return(reticulate::py_to_r(np$array(mat)))
        }
        return(NULL)
    }
    
    se_py <- get_obsm_matrix("se")
    zscore_py <- get_obsm_matrix("zscore")
    pvalue_py <- get_obsm_matrix("pvalue")

    ## ----------------------------------------------------------------
    ## 3. Transpose to R format (proteins x cells)
    ## ----------------------------------------------------------------
    beta <- t(beta_py)
    rownames(beta) <- protein_names
    colnames(beta) <- cell_names
    
    if (!is.null(se_py)) {
        se <- t(se_py)
        rownames(se) <- protein_names
        colnames(se) <- cell_names
    } else {
        se <- NULL
    }
    
    if (!is.null(zscore_py)) {
        zscore <- t(zscore_py)
        rownames(zscore) <- protein_names
        colnames(zscore) <- cell_names
    } else {
        zscore <- NULL
    }
    
    if (!is.null(pvalue_py)) {
        pvalue <- t(pvalue_py)
        rownames(pvalue) <- protein_names
        colnames(pvalue) <- cell_names
    } else {
        pvalue <- NULL
    }
    
    cat("Output shape (R format):", nrow(beta), "x", ncol(beta), "(proteins x cells)\n")

    ## ----------------------------------------------------------------
    ## 4. Create Seurat object
    ## ----------------------------------------------------------------
    old_option <- getOption("Seurat.object.assay.version")
    options(Seurat.object.assay.version = "v3")
    on.exit(options(Seurat.object.assay.version = old_option), add = TRUE)
    
    seu <- Seurat::CreateSeuratObject(
        counts = beta,
        assay = "SecAct"
    )
    
    ## Store SecAct results in misc
    seu@misc$SecAct_output$SecretedProteinActivity <- list(
        beta   = beta,
        se     = se,
        zscore = zscore,
        pvalue = pvalue
    )
    seu@misc$protein_names <- protein_names

    cat("Created Seurat object with SecAct results\n")
    
    return(seu)
}
