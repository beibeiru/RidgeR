#' =============================================================================
#' RidgeR: Secreted Protein Activity Inference
#' =============================================================================
#'
#' This module provides R wrappers for ridge regression-based activity inference.
#' Multiple implementations are available with different performance characteristics:
#'
#'   - Legacy:  Original .C interface (single-threaded)
#'   - Old:     .Call interface with Y-permutation (single-threaded)
#'   - Old2:    .Call interface with T-permutation (single-threaded, cache-friendly)
#'   - New:     .Call interface with Y-permutation (multi-threaded, OpenMP)
#'   - New2:    .Call interface with T-permutation (multi-threaded, OpenMP)
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
#' @return List with aligned and scaled X and Y matrices
#' @keywords internal
.prepare_matrices <- function(X, Y) {
    # Find overlapping genes
    common_genes <- intersect(rownames(Y), rownames(X))

    if (length(common_genes) == 0) {
        stop("No common genes found between signature and expression matrices")
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
#' @param v Flat result vector
#' @param X Signature matrix (for column names = proteins)
#' @param Y Expression matrix (for column names = samples)
#' @return Matrix with appropriate dimensions and names
#' @keywords internal
.format_result <- function(v, X, Y) {
    matrix(
        v,
        nrow = ncol(X),
        ncol = ncol(Y),
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

    # Call C function
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

    # Format and return results
    list(
        beta   = .format_result(res$beta,   X, Y),
        se     = .format_result(res$se,     X, Y),
        zscore = .format_result(res$zscore, X, Y),
        pvalue = .format_result(res$pvalue, X, Y)
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
SecAct.inference.gsl.old <- function(Y,
                                      SigMat = "SecAct",
                                      lambda = 5e+05,
                                      nrand = 1000) {
    # Load and prepare data
    X <- .load_signature(SigMat)
    prep <- .prepare_matrices(X, Y)

    # Call C function
    res <- .Call(
        "ridgeReg_old_interface",
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
SecAct.inference.gsl.old2 <- function(Y,
                                       SigMat = "SecAct",
                                       lambda = 5e+05,
                                       nrand = 1000) {
    # Load and prepare data
    X <- .load_signature(SigMat)
    prep <- .prepare_matrices(X, Y)

    # Call C function
    res <- .Call(
        "ridgeRegTperm_old_interface",
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
SecAct.inference.gsl.new <- function(Y,
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
        "ridgeRegFast_interface",
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
#' @inheritParams SecAct.inference.gsl.new
#'
#' @return Named list with beta, se, zscore, pvalue matrices (proteins x samples)
#'
#' @export
SecAct.inference.gsl.new2 <- function(Y,
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
        "ridgeRegTperm_interface",
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
# SECTION 5: TESTING & COMPARISON UTILITIES
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
#' @param method Which implementation to use: "legacy", "old", "old2", "new", "new2"
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
    method = "new"
) {
    # Compute differential expression if needed
    if (is.differential) {
        Y <- inputProfile
        if (ncol(Y) == 1) colnames(Y) <- "Change"
    } else {
        if (is.null(inputProfile_control)) {
            Y <- inputProfile - rowMeans(inputProfile)
        } else {
            if (is.paired) {
                Y <- inputProfile - inputProfile_control[, colnames(inputProfile), drop = FALSE]
            } else {
                Y <- inputProfile - rowMeans(inputProfile_control)
            }
            
            if (!is.singleSampleLevel) {
                Y <- matrix(rowMeans(Y), ncol = 1, 
                           dimnames = list(rownames(Y), "Change"))
            }
        }
    }
    
    # Load signature matrix
    X <- .load_signature(sigMatrix)
    
    # Filter signatures if requested
    if (sigFilter) {
        X <- X[, colnames(X) %in% rownames(Y), drop = FALSE]
    }
    
    # Group similar signatures if requested
    if (is.group.sig) {
        X <- .group_signatures(X, is.group.cor)
    }
    
    # Find overlapping genes
    olp <- intersect(rownames(Y), rownames(X))
    
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
    
    # Run inference using selected method
    ncores <- .default_ncores(ncores)
    
    res <- switch(method,
        "legacy" = {
            n <- nrow(X)
            p <- ncol(X)
            m <- ncol(Y)
            len <- p * m
            
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
        "old" = {
            raw <- .Call(
                "ridgeReg_old_interface",
                t(X), t(Y),
                as.numeric(lambda),
                as.integer(nrand),
                PACKAGE = "RidgeR"
            )
            list(
                beta   = matrix(raw$beta, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        "old2" = {
            raw <- .Call(
                "ridgeRegTperm_old_interface",
                t(X), t(Y),
                as.numeric(lambda),
                as.integer(nrand),
                PACKAGE = "RidgeR"
            )
            list(
                beta   = matrix(raw$beta, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        "new" = {
            raw <- .Call(
                "ridgeRegFast_interface",
                t(X), t(Y),
                as.numeric(lambda),
                as.integer(nrand),
                as.integer(ncores),
                PACKAGE = "RidgeR"
            )
            list(
                beta   = matrix(raw$beta, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        "new2" = {
            raw <- .Call(
                "ridgeRegTperm_interface",
                t(X), t(Y),
                as.numeric(lambda),
                as.integer(nrand),
                as.integer(ncores),
                PACKAGE = "RidgeR"
            )
            list(
                beta   = matrix(raw$beta, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                se     = matrix(raw$se, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                zscore = matrix(raw$zscore, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y))),
                pvalue = matrix(raw$pvalue, nrow = ncol(X), ncol = ncol(Y),
                               dimnames = list(colnames(X), colnames(Y)))
            )
        },
        stop("Unknown method: ", method)
    )
    
    # Expand grouped signatures back to individual rows
    if (is.group.sig) {
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
    
    res
}


#' Compare All Implementation Methods
#'
#' Utility function to verify that all implementations produce identical
#' (within tolerance) results and to compare their execution times.
#'
#' @inheritParams SecAct.inference.gsl.new
#' @param tol Tolerance for numerical comparison (default: 1e-10)
#'
#' @return Invisibly returns a named list containing results from all methods:
#'   \itemize{
#'     \item \code{legacy}: Results from .C legacy implementation
#'     \item \code{old}: Results from single-threaded Y-permutation
#'     \item \code{old2}: Results from single-threaded T-permutation
#'     \item \code{new}: Results from multi-threaded Y-permutation
#'     \item \code{new2}: Results from multi-threaded T-permutation
#'   }
#'
#' @examples
#' \dontrun{
#' # Create test data
#' set.seed(123)
#' Y <- matrix(rnorm(1000 * 50), nrow = 1000, ncol = 50)
#' rownames(Y) <- paste0("Gene", 1:1000)
#' colnames(Y) <- paste0("Sample", 1:50)
#'
#' # Compare methods with small nrand for quick testing
#' results <- SecAct.compare.methods(Y, nrand = 100)
#' }
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

    cat("Running gsl.old (single-threaded, Y-permutation)...\n")
    t_old <- time_it({
        res_old <- SecAct.inference.gsl.old(Y, SigMat, lambda, nrand)
    })

    cat("Running gsl.old2 (single-threaded, T-permutation)...\n")
    t_old2 <- time_it({
        res_old2 <- SecAct.inference.gsl.old2(Y, SigMat, lambda, nrand)
    })

    cat("Running gsl.new (multi-threaded, Y-permutation)...\n")
    t_new <- time_it({
        res_new <- SecAct.inference.gsl.new(Y, SigMat, lambda, nrand, ncores)
    })

    cat("Running gsl.new2 (multi-threaded, T-permutation)...\n")
    t_new2 <- time_it({
        res_new2 <- SecAct.inference.gsl.new2(Y, SigMat, lambda, nrand, ncores)
    })

    # Helper to check equality within tolerance
    check_equal <- function(a, b) {
        max_diff <- max(abs(a - b), na.rm = TRUE)
        max_diff < tol
    }

    # Report consistency
    cat("\n")
    cat("=== Consistency Check (z-score comparison) ===\n")
    cat(sprintf("  legacy vs old  : %s\n",
                ifelse(check_equal(res_legacy$zscore, res_old$zscore), "PASS", "FAIL")))
    cat(sprintf("  old    vs new  : %s\n",
                ifelse(check_equal(res_old$zscore, res_new$zscore), "PASS", "FAIL")))
    cat(sprintf("  old2   vs new2 : %s\n",
                ifelse(check_equal(res_old2$zscore, res_new2$zscore), "PASS", "FAIL")))
    cat(sprintf("  new    vs new2 : %s\n",
                ifelse(check_equal(res_new$zscore, res_new2$zscore), "PASS", "FAIL")))

    # Report timing
    cat("\n")
    cat("=== Timing Summary ===\n")
    cat(sprintf("  gsl.legacy : %6.2f seconds\n", t_legacy))
    cat(sprintf("  gsl.old    : %6.2f seconds\n", t_old))
    cat(sprintf("  gsl.old2   : %6.2f seconds\n", t_old2))
    cat(sprintf("  gsl.new    : %6.2f seconds\n", t_new))
    cat(sprintf("  gsl.new2   : %6.2f seconds\n", t_new2))

    # Return all results invisibly
    invisible(list(
        legacy = res_legacy,
        old    = res_old,
        old2   = res_old2,
        new    = res_new,
        new2   = res_new2
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
    
    cat("Method: old (single-threaded, Y-permutation)...\n")
    t_old <- time_it({
        res_old <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, method = "old"
        )
    })
    
    cat("Method: old2 (single-threaded, T-permutation)...\n")
    t_old2 <- time_it({
        res_old2 <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, method = "old2"
        )
    })
    
    cat("Method: new (multi-threaded, Y-permutation)...\n")
    t_new <- time_it({
        res_new <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, ncores = ncores, method = "new"
        )
    })
    
    cat("Method: new2 (multi-threaded, T-permutation)...\n")
    t_new2 <- time_it({
        res_new2 <- SecAct.activity.inference(
            inputProfile, is.differential = is.differential,
            sigMatrix = sigMatrix, is.group.sig = is.group.sig,
            is.group.cor = is.group.cor, lambda = lambda,
            nrand = nrand, ncores = ncores, method = "new2"
        )
    })
    
    # Report results
    cat("\n=== Z-score Comparison (first 6 rows) ===\n")
    cat("\nLegacy:\n")
    print(head(res_legacy$zscore))
    cat("\nNew:\n")
    print(head(res_new$zscore))
    
    # Check consistency
    tol <- 1e-10
    check_equal <- function(a, b) {
        max(abs(a - b), na.rm = TRUE) < tol
    }
    
    cat("\n=== Consistency Check ===\n")
    cat(sprintf("  legacy vs old  : %s\n",
                ifelse(check_equal(res_legacy$zscore, res_old$zscore), "PASS", "FAIL")))
    cat(sprintf("  old    vs new  : %s\n",
                ifelse(check_equal(res_old$zscore, res_new$zscore), "PASS", "FAIL")))
    cat(sprintf("  old2   vs new2 : %s\n",
                ifelse(check_equal(res_old2$zscore, res_new2$zscore), "PASS", "FAIL")))
    cat(sprintf("  new    vs new2 : %s\n",
                ifelse(check_equal(res_new$zscore, res_new2$zscore), "PASS", "FAIL")))
    
    cat("\n=== Timing Summary ===\n")
    cat(sprintf("  legacy : %6.2f seconds\n", t_legacy))
    cat(sprintf("  old    : %6.2f seconds\n", t_old))
    cat(sprintf("  old2   : %6.2f seconds\n", t_old2))
    cat(sprintf("  new    : %6.2f seconds\n", t_new))
    cat(sprintf("  new2   : %6.2f seconds\n", t_new2))
    
    invisible(list(
        legacy = res_legacy,
        old    = res_old,
        old2   = res_old2,
        new    = res_new,
        new2   = res_new2
    ))
}
