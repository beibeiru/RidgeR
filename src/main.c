/**
 * =============================================================================
 * RidgeR: Ridge Regression with Permutation Testing
 * =============================================================================
 *
 * This module implements ridge regression with permutation-based significance
 * testing for secreted protein activity inference. Multiple implementations
 * are provided for different use cases:
 *
 *   Version 0 (Legacy):  .C interface, single-threaded, Y-permutation
 *   Version 1 (Old):     .Call interface, single-threaded, Y-permutation
 *   Version 2 (Old2):    .Call interface, single-threaded, T-permutation
 *   Version 3 (New):     .Call interface, multi-threaded (OpenMP), Y-permutation
 *   Version 4 (New2):    .Call interface, multi-threaded (OpenMP), T-permutation
 *
 * Dependencies: GSL (GNU Scientific Library), OpenMP (optional)
 *
 * =============================================================================
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <string.h>
#include <math.h>

#ifdef _OPENMP
    #include <omp.h>
#endif


/* =============================================================================
 * SECTION 1: UTILITY FUNCTIONS
 * =============================================================================
 */

/**
 * Fisher-Yates shuffle algorithm (in-place)
 *
 * @param array  Integer array to shuffle
 * @param n      Length of array
 */
static void shuffle_array(int *array, const int n) {
    for (int i = 0; i < n - 1; i++) {
        int j = i + rand() / (RAND_MAX / (n - i) + 1);
        int tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}


/**
 * Generate a table of random permutations
 *
 * Creates nrand independent permutations of indices [0, n-1].
 * Uses a fixed seed (0) for reproducibility.
 *
 * @param n      Number of elements per permutation
 * @param nrand  Number of permutations to generate
 * @return       Flat array of size nrand * n (caller must free)
 */
static int* generate_permutation_table(int n, int nrand) {
    size_t table_size = (size_t)nrand * (size_t)n;
    int *table = (int*)malloc(table_size * sizeof(int));
    int *base_indices = (int*)malloc((size_t)n * sizeof(int));

    /* Initialize base index array */
    for (int i = 0; i < n; i++) {
        base_indices[i] = i;
    }

    /* Generate permutations with fixed seed for reproducibility */
    srand(0);
    for (int perm = 0; perm < nrand; perm++) {
        shuffle_array(base_indices, n);
        memcpy(table + ((size_t)perm * n), base_indices, n * sizeof(int));
    }

    free(base_indices);
    return table;
}


/* =============================================================================
 * SECTION 2: GSL MATRIX HELPERS
 * =============================================================================
 */

/**
 * Create a GSL matrix view wrapping an R vector (no data copy)
 *
 * Note: The returned matrix takes ownership flag but does NOT own the data.
 *       Use gsl_matrix_partial_free() to free only the wrapper.
 *
 * @param vec  Pointer to column-major R data
 * @param nr   Number of rows
 * @param nc   Number of columns
 * @return     GSL matrix wrapper (caller must call gsl_matrix_partial_free)
 */
static gsl_matrix* wrap_r_vector_as_gsl_matrix(double *vec, size_t nr, size_t nc) {
    gsl_block *block = (gsl_block*)malloc(sizeof(gsl_block));
    gsl_matrix *mat = (gsl_matrix*)malloc(sizeof(gsl_matrix));

    mat->size1 = nr;
    mat->size2 = nc;
    mat->tda = nc;
    mat->owner = 1;  /* Indicates we own the block struct, not the data */

    block->data = vec;
    block->size = nr * nc;

    mat->data = vec;
    mat->block = block;

    return mat;
}


/**
 * Free GSL matrix wrapper without freeing the underlying data
 *
 * Use this for matrices created with wrap_r_vector_as_gsl_matrix()
 */
static void gsl_matrix_partial_free(gsl_matrix *mat) {
    if (mat) {
        if (mat->block) {
            free(mat->block);
        }
        free(mat);
    }
}


/* =============================================================================
 * SECTION 3: R INTERFACE HELPERS
 * =============================================================================
 */

/**
 * Create the standard result list: (beta, se, zscore, pvalue)
 *
 * @param beta    SEXP for beta coefficients
 * @param se      SEXP for standard errors
 * @param zscore  SEXP for z-scores
 * @param pvalue  SEXP for p-values
 * @return        Named list with 4 elements
 */
static SEXP create_result_list(SEXP beta, SEXP se, SEXP zscore, SEXP pvalue) {
    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SEXP names = PROTECT(allocVector(STRSXP, 4));

    SET_VECTOR_ELT(result, 0, beta);
    SET_VECTOR_ELT(result, 1, se);
    SET_VECTOR_ELT(result, 2, zscore);
    SET_VECTOR_ELT(result, 3, pvalue);

    SET_STRING_ELT(names, 0, mkChar("beta"));
    SET_STRING_ELT(names, 1, mkChar("se"));
    SET_STRING_ELT(names, 2, mkChar("zscore"));
    SET_STRING_ELT(names, 3, mkChar("pvalue"));

    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(2);
    return result;
}


/**
 * Compute the ridge regression projection matrix T = (X'X + λI)^{-1} X'
 *
 * @param Xt      Transposed design matrix X' (p x n)
 * @param lambda  Ridge penalty parameter
 * @param I_out   Output: p x p identity matrix (modified in place)
 * @param T_out   Output: Projection matrix T (p x n)
 */
static void compute_projection_matrix(const gsl_matrix *Xt, double lambda,
                                       gsl_matrix *I_out, gsl_matrix *T_out) {
    size_t p = Xt->size1;

    /* I = λI + X'X */
    gsl_matrix_set_identity(I_out);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_out);

    /* I = (X'X + λI)^{-1} via Cholesky decomposition */
    gsl_linalg_cholesky_decomp(I_out);
    gsl_linalg_cholesky_invert(I_out);

    /* T = (X'X + λI)^{-1} X' */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_out, Xt, 0.0, T_out);
}


/**
 * Finalize permutation statistics
 *
 * Converts accumulated sums to standard errors, z-scores, and p-values.
 *
 * @param bv      Beta coefficients (observed)
 * @param sv      Sum of squares (input) / Standard error (output)
 * @param zv      Sum (input) / Z-score (output)
 * @param pv      Count of |perm| >= |obs| (input) / P-value (output)
 * @param len     Number of elements
 * @param nrand   Number of permutations
 */
static void finalize_permutation_stats(const double *bv, double *sv, double *zv,
                                        double *pv, R_xlen_t len, int nrand) {
    double inv_nrand = 1.0 / nrand;

    for (R_xlen_t i = 0; i < len; i++) {
        double mean = zv[i] * inv_nrand;
        double var = (sv[i] * inv_nrand) - (mean * mean);
        double std = sqrt(fmax(0.0, var));

        sv[i] = std;
        zv[i] = (std > 1e-12) ? (bv[i] - mean) / std : 0.0;
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
    }
}


/* =============================================================================
 * SECTION 4: VERSION 0 - LEGACY .C INTERFACE
 * =============================================================================
 *
 * Original implementation using R's .C() interface.
 * Single-threaded with Y-row permutation strategy.
 */

void ridgeReg(double *X_vec, double *Y_vec,
              int *n_ptr, int *p_ptr, int *m_ptr,
              double *lambda_ptr, double *nrand_ptr,
              double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec) {

    /* Extract dimensions */
    size_t n = (size_t)*n_ptr;  /* Number of samples */
    size_t p = (size_t)*p_ptr;  /* Number of predictors (proteins) */
    size_t m = (size_t)*m_ptr;  /* Number of responses (samples in Y) */
    int nrand = (int)*nrand_ptr;
    double lambda = *lambda_ptr;
    size_t total_elements = p * m;

    /* Wrap R vectors as GSL matrices */
    gsl_matrix *Xt = wrap_r_vector_as_gsl_matrix(X_vec, p, n);   /* X transposed */
    gsl_matrix *Yt = wrap_r_vector_as_gsl_matrix(Y_vec, m, n);   /* Y transposed */
    gsl_matrix *beta = wrap_r_vector_as_gsl_matrix(beta_vec, p, m);

    /* Allocate working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Yr = gsl_matrix_alloc(n, m);      /* Y in row-major form */
    gsl_matrix *Yp = gsl_matrix_alloc(n, m);      /* Permuted Y */
    gsl_matrix *br = gsl_matrix_alloc(p, m);      /* Permuted beta */
    gsl_matrix *sum_beta = gsl_matrix_alloc(p, m);
    gsl_matrix *sum_beta_sq = gsl_matrix_alloc(p, m);

    /* Compute projection matrix and observed beta */
    compute_projection_matrix(Xt, lambda, I, T);
    gsl_matrix_transpose_memcpy(Yr, Yt);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, beta);

    /* Initialize accumulators */
    int *perm_table = generate_permutation_table((int)n, nrand);
    gsl_matrix_set_zero(sum_beta);
    gsl_matrix_set_zero(sum_beta_sq);
    for (size_t i = 0; i < total_elements; i++) {
        pvalue_vec[i] = 0.0;
    }

    /* Permutation loop */
    for (int perm = 0; perm < nrand; perm++) {
        int *perm_indices = perm_table + (perm * n);

        /* Apply permutation to Y rows */
        for (size_t row = 0; row < n; row++) {
            memcpy(Yp->data + (row * m),
                   Yr->data + ((size_t)perm_indices[row] * m),
                   m * sizeof(double));
        }

        /* Compute permuted coefficients */
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yp, 0.0, br);

        /* Accumulate statistics */
        for (size_t j = 0; j < total_elements; j++) {
            if (fabs(br->data[j]) >= fabs(beta->data[j])) {
                pvalue_vec[j]++;
            }
        }
        gsl_matrix_add(sum_beta, br);
        gsl_matrix_mul_elements(br, br);
        gsl_matrix_add(sum_beta_sq, br);
    }

    /* Finalize statistics */
    double inv_nrand = 1.0 / nrand;
    for (size_t i = 0; i < total_elements; i++) {
        pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
        double mean = sum_beta->data[i] * inv_nrand;
        double var = (sum_beta_sq->data[i] * inv_nrand) - (mean * mean);
        double std = sqrt(fmax(0.0, var));
        se_vec[i] = std;
        zscore_vec[i] = (std > 1e-12) ? (beta_vec[i] - mean) / std : 0.0;
    }

    /* Cleanup */
    free(perm_table);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Yr);
    gsl_matrix_free(Yp);
    gsl_matrix_free(br);
    gsl_matrix_free(sum_beta);
    gsl_matrix_free(sum_beta_sq);
    gsl_matrix_partial_free(Xt);
    gsl_matrix_partial_free(Yt);
    gsl_matrix_partial_free(beta);
}


/* =============================================================================
 * SECTION 5: VERSION 1 - SINGLE-THREADED Y-PERMUTATION (.Call)
 * =============================================================================
 */

SEXP ridgeReg_old_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp) {

    /* Extract dimensions */
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0];
    size_t p = (size_t)INTEGER(x_dim)[1];
    size_t m = (size_t)INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];
    int nrand = asInteger(nrand_sexp);
    double lambda = asReal(lambda_sexp);
    R_xlen_t len = (R_xlen_t)p * m;

    /* Allocate output vectors */
    SEXP beta_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP se_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP zscore_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP pvalue_sexp = PROTECT(allocVector(REALSXP, len));

    double *bv = REAL(beta_sexp);
    double *sv = REAL(se_sexp);
    double *zv = REAL(zscore_sexp);
    double *pv = REAL(pvalue_sexp);

    memset(pv, 0, len * sizeof(double));
    memset(sv, 0, len * sizeof(double));
    memset(zv, 0, len * sizeof(double));

    /* Create GSL views of input matrices */
    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);

    /* Allocate working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Yr = gsl_matrix_alloc(n, m);
    gsl_matrix *Yp = gsl_matrix_alloc(n, m);
    gsl_matrix *br = gsl_matrix_alloc(p, m);

    /* Compute projection matrix and observed beta */
    compute_projection_matrix(&Xt.matrix, lambda, I, T);
    gsl_matrix_transpose_memcpy(Yr, &Yt.matrix);
    gsl_matrix_view beta_view = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, &beta_view.matrix);

    /* Permutation testing */
    int *perm_table = generate_permutation_table((int)n, nrand);

    for (int perm = 0; perm < nrand; perm++) {
        int *perm_indices = perm_table + (perm * n);

        /* Apply permutation */
        for (size_t row = 0; row < n; row++) {
            memcpy(Yp->data + (row * m),
                   Yr->data + ((size_t)perm_indices[row] * m),
                   m * sizeof(double));
        }

        /* Compute permuted coefficients */
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yp, 0.0, br);

        /* Accumulate statistics */
        for (R_xlen_t k = 0; k < len; k++) {
            if (fabs(br->data[k]) >= fabs(bv[k])) pv[k]++;
            zv[k] += br->data[k];
            sv[k] += br->data[k] * br->data[k];
        }
    }

    /* Finalize statistics */
    finalize_permutation_stats(bv, sv, zv, pv, len, nrand);

    /* Cleanup */
    free(perm_table);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Yr);
    gsl_matrix_free(Yp);
    gsl_matrix_free(br);

    SEXP result = create_result_list(beta_sexp, se_sexp, zscore_sexp, pvalue_sexp);
    UNPROTECT(4);
    return result;
}


/* =============================================================================
 * SECTION 6: VERSION 2 - SINGLE-THREADED T-PERMUTATION (.Call)
 * =============================================================================
 *
 * Uses T-column permutation (scatter approach) which can be more cache-friendly
 * for certain matrix sizes.
 */

SEXP ridgeRegTperm_old_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp) {

    /* Extract dimensions */
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0];
    size_t p = (size_t)INTEGER(x_dim)[1];
    size_t m = (size_t)INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];
    int nrand = asInteger(nrand_sexp);
    double lambda = asReal(lambda_sexp);
    R_xlen_t len = (R_xlen_t)p * m;

    /* Allocate output vectors */
    SEXP beta_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP se_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP zscore_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP pvalue_sexp = PROTECT(allocVector(REALSXP, len));

    double *bv = REAL(beta_sexp);
    double *sv = REAL(se_sexp);
    double *zv = REAL(zscore_sexp);
    double *pv = REAL(pvalue_sexp);

    memset(pv, 0, len * sizeof(double));
    memset(sv, 0, len * sizeof(double));
    memset(zv, 0, len * sizeof(double));

    /* Create GSL views */
    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);

    /* Allocate working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Yr = gsl_matrix_alloc(n, m);
    gsl_matrix *Tt_original = gsl_matrix_alloc(n, p);  /* T transposed */
    gsl_matrix *Tt_permuted = gsl_matrix_alloc(n, p);
    gsl_matrix *br = gsl_matrix_alloc(p, m);

    /* Compute projection matrix and observed beta */
    compute_projection_matrix(&Xt.matrix, lambda, I, T);
    gsl_matrix_transpose_memcpy(Yr, &Yt.matrix);
    gsl_matrix_view beta_view = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, &beta_view.matrix);

    /* Prepare T transpose for permutation */
    gsl_matrix_transpose_memcpy(Tt_original, T);

    /* Permutation testing (permute T columns via scatter) */
    int *perm_table = generate_permutation_table((int)n, nrand);

    for (int perm = 0; perm < nrand; perm++) {
        int *perm_indices = perm_table + (perm * n);

        /* Scatter permutation: Tt_permuted[perm[i], :] = Tt_original[i, :] */
        for (size_t i = 0; i < n; i++) {
            memcpy(Tt_permuted->data + ((size_t)perm_indices[i] * p),
                   Tt_original->data + (i * p),
                   p * sizeof(double));
        }

        /* Compute permuted coefficients: beta = Tt_permuted' * Yr */
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Tt_permuted, Yr, 0.0, br);

        /* Accumulate statistics */
        for (R_xlen_t k = 0; k < len; k++) {
            if (fabs(br->data[k]) >= fabs(bv[k])) pv[k]++;
            zv[k] += br->data[k];
            sv[k] += br->data[k] * br->data[k];
        }
    }

    /* Finalize statistics */
    finalize_permutation_stats(bv, sv, zv, pv, len, nrand);

    /* Cleanup */
    free(perm_table);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Yr);
    gsl_matrix_free(Tt_original);
    gsl_matrix_free(Tt_permuted);
    gsl_matrix_free(br);

    SEXP result = create_result_list(beta_sexp, se_sexp, zscore_sexp, pvalue_sexp);
    UNPROTECT(4);
    return result;
}


/* =============================================================================
 * SECTION 7: VERSION 3 - MULTI-THREADED Y-PERMUTATION (OpenMP)
 * =============================================================================
 */

SEXP ridgeRegFast_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp,
                            SEXP nrand_sexp, SEXP ncores_sexp) {

    /* Configure OpenMP threads */
    int ncores = asInteger(ncores_sexp);
    #ifdef _OPENMP
        if (ncores > 0) omp_set_num_threads(ncores);
    #endif

    /* Extract dimensions */
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0];
    size_t p = (size_t)INTEGER(x_dim)[1];
    size_t m = (size_t)INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];
    int nrand = asInteger(nrand_sexp);
    double lambda = asReal(lambda_sexp);
    R_xlen_t len = (R_xlen_t)p * m;

    /* Allocate output vectors */
    SEXP beta_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP se_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP zscore_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP pvalue_sexp = PROTECT(allocVector(REALSXP, len));

    double *bv = REAL(beta_sexp);
    double *sv = REAL(se_sexp);
    double *zv = REAL(zscore_sexp);
    double *pv = REAL(pvalue_sexp);

    memset(pv, 0, len * sizeof(double));
    memset(sv, 0, len * sizeof(double));
    memset(zv, 0, len * sizeof(double));

    /* Create GSL views */
    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);

    /* Allocate shared working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Yr = gsl_matrix_alloc(n, m);

    /* Compute projection matrix and observed beta */
    compute_projection_matrix(&Xt.matrix, lambda, I, T);
    gsl_matrix_transpose_memcpy(Yr, &Yt.matrix);
    gsl_matrix_view beta_view = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, &beta_view.matrix);

    /* Prepare permutation table */
    int *perm_table = generate_permutation_table((int)n, nrand);

    /* Determine actual thread count */
    int actual_threads = 1;
    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        actual_threads = omp_get_num_threads();
    }
    #endif

    /* Allocate thread-local accumulators */
    double *thread_sums = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *thread_sumsq = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *thread_counts = (double*)calloc((size_t)actual_threads * len, sizeof(double));

    /* Parallel permutation loop */
    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
            tid = omp_get_thread_num();
        #endif

        /* Thread-local accumulator pointers */
        double *my_sums = &thread_sums[(size_t)tid * len];
        double *my_sumsq = &thread_sumsq[(size_t)tid * len];
        double *my_counts = &thread_counts[(size_t)tid * len];

        /* Thread-local working matrices */
        gsl_matrix *Yp = gsl_matrix_alloc(n, m);
        gsl_matrix *br = gsl_matrix_alloc(p, m);

        #pragma omp for schedule(dynamic)
        for (int perm = 0; perm < nrand; perm++) {
            int *perm_indices = perm_table + ((size_t)perm * n);

            /* Apply permutation */
            for (size_t row = 0; row < n; row++) {
                memcpy(Yp->data + (row * m),
                       Yr->data + ((size_t)perm_indices[row] * m),
                       m * sizeof(double));
            }

            /* Compute permuted coefficients */
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yp, 0.0, br);

            /* Accumulate statistics */
            for (R_xlen_t k = 0; k < len; k++) {
                if (fabs(br->data[k]) >= fabs(bv[k])) my_counts[k]++;
                my_sums[k] += br->data[k];
                my_sumsq[k] += br->data[k] * br->data[k];
            }
        }

        gsl_matrix_free(Yp);
        gsl_matrix_free(br);
    }

    /* Reduce thread-local results */
    for (int t = 0; t < actual_threads; t++) {
        size_t offset = (size_t)t * len;
        for (R_xlen_t i = 0; i < len; i++) {
            zv[i] += thread_sums[offset + i];
            sv[i] += thread_sumsq[offset + i];
            pv[i] += thread_counts[offset + i];
        }
    }

    /* Finalize statistics */
    finalize_permutation_stats(bv, sv, zv, pv, len, nrand);

    /* Cleanup */
    free(thread_sums);
    free(thread_sumsq);
    free(thread_counts);
    free(perm_table);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Yr);

    SEXP result = create_result_list(beta_sexp, se_sexp, zscore_sexp, pvalue_sexp);
    UNPROTECT(4);
    return result;
}


/* =============================================================================
 * SECTION 8: VERSION 4 - MULTI-THREADED T-PERMUTATION (OpenMP)
 * =============================================================================
 */

SEXP ridgeRegTperm_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp,
                              SEXP nrand_sexp, SEXP ncores_sexp) {

    /* Configure OpenMP threads */
    int ncores = asInteger(ncores_sexp);
    #ifdef _OPENMP
        if (ncores > 0) omp_set_num_threads(ncores);
    #endif

    /* Extract dimensions */
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0];
    size_t p = (size_t)INTEGER(x_dim)[1];
    size_t m = (size_t)INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];
    int nrand = asInteger(nrand_sexp);
    double lambda = asReal(lambda_sexp);
    R_xlen_t len = (R_xlen_t)p * m;

    /* Allocate output vectors */
    SEXP beta_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP se_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP zscore_sexp = PROTECT(allocVector(REALSXP, len));
    SEXP pvalue_sexp = PROTECT(allocVector(REALSXP, len));

    double *bv = REAL(beta_sexp);
    double *sv = REAL(se_sexp);
    double *zv = REAL(zscore_sexp);
    double *pv = REAL(pvalue_sexp);

    memset(pv, 0, len * sizeof(double));
    memset(sv, 0, len * sizeof(double));
    memset(zv, 0, len * sizeof(double));

    /* Create GSL views */
    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);

    /* Allocate shared working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Yr = gsl_matrix_alloc(n, m);
    gsl_matrix *Tt_original = gsl_matrix_alloc(n, p);

    /* Compute projection matrix and observed beta */
    compute_projection_matrix(&Xt.matrix, lambda, I, T);
    gsl_matrix_transpose_memcpy(Yr, &Yt.matrix);
    gsl_matrix_view beta_view = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, &beta_view.matrix);

    /* Prepare T transpose */
    gsl_matrix_transpose_memcpy(Tt_original, T);

    /* Prepare permutation table */
    int *perm_table = generate_permutation_table((int)n, nrand);

    /* Determine actual thread count */
    int actual_threads = 1;
    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        actual_threads = omp_get_num_threads();
    }
    #endif

    /* Allocate thread-local accumulators */
    double *thread_sums = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *thread_sumsq = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *thread_counts = (double*)calloc((size_t)actual_threads * len, sizeof(double));

    /* Parallel permutation loop */
    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
            tid = omp_get_thread_num();
        #endif

        /* Thread-local accumulator pointers */
        double *my_sums = &thread_sums[(size_t)tid * len];
        double *my_sumsq = &thread_sumsq[(size_t)tid * len];
        double *my_counts = &thread_counts[(size_t)tid * len];

        /* Thread-local working matrices */
        gsl_matrix *Tt_permuted = gsl_matrix_alloc(n, p);
        gsl_matrix *br = gsl_matrix_alloc(p, m);

        #pragma omp for schedule(dynamic, 16)
        for (int perm = 0; perm < nrand; perm++) {
            int *perm_indices = perm_table + ((size_t)perm * n);

            /* Scatter permutation */
            for (size_t i = 0; i < n; i++) {
                memcpy(Tt_permuted->data + ((size_t)perm_indices[i] * p),
                       Tt_original->data + (i * p),
                       p * sizeof(double));
            }

            /* Compute permuted coefficients */
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Tt_permuted, Yr, 0.0, br);

            /* Accumulate statistics */
            for (R_xlen_t k = 0; k < len; k++) {
                if (fabs(br->data[k]) >= fabs(bv[k])) my_counts[k]++;
                my_sums[k] += br->data[k];
                my_sumsq[k] += br->data[k] * br->data[k];
            }
        }

        gsl_matrix_free(Tt_permuted);
        gsl_matrix_free(br);
    }

    /* Reduce thread-local results */
    for (int t = 0; t < actual_threads; t++) {
        size_t offset = (size_t)t * len;
        for (R_xlen_t i = 0; i < len; i++) {
            zv[i] += thread_sums[offset + i];
            sv[i] += thread_sumsq[offset + i];
            pv[i] += thread_counts[offset + i];
        }
    }

    /* Finalize statistics */
    finalize_permutation_stats(bv, sv, zv, pv, len, nrand);

    /* Cleanup */
    free(thread_sums);
    free(thread_sumsq);
    free(thread_counts);
    free(perm_table);
    gsl_matrix_free(Tt_original);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Yr);

    SEXP result = create_result_list(beta_sexp, se_sexp, zscore_sexp, pvalue_sexp);
    UNPROTECT(4);
    return result;
}


/* =============================================================================
 * SECTION 9: R REGISTRATION
 * =============================================================================
 */

static const R_CMethodDef cMethods[] = {
    {"ridgeReg", (DL_FUNC) &ridgeReg, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef callMethods[] = {
    {"ridgeReg_old_interface",      (DL_FUNC) &ridgeReg_old_interface,      4},
    {"ridgeRegTperm_old_interface", (DL_FUNC) &ridgeRegTperm_old_interface, 4},
    {"ridgeRegFast_interface",      (DL_FUNC) &ridgeRegFast_interface,      5},
    {"ridgeRegTperm_interface",     (DL_FUNC) &ridgeRegTperm_interface,     5},
    {NULL, NULL, 0}
};

void R_init_RidgeR(DllInfo *dll) {
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
