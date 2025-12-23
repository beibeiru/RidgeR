/**
 * =============================================================================
 * RidgeR: Ridge Regression with Permutation Testing
 * =============================================================================
 *
 * This module implements ridge regression with permutation-based significance
 * testing for secreted protein activity inference. Multiple implementations
 * are provided for different use cases:
 *
 *   legacy: .C interface, single-threaded, Y-permutation
 *   styp:   .Call interface, Single-Thread Y-Permutation
 *   sttp:   .Call interface, Single-Thread T-Permutation (cache-friendly)
 *   mtyp:   .Call interface, Multi-Thread Y-Permutation (OpenMP)
 *   mttp:   .Call interface, Multi-Thread T-Permutation (OpenMP, cache-friendly)
 *
 * Dependencies: GSL (GNU Scientific Library), OpenMP (optional)
 *
 * RNG: All versions use GSL's Mersenne Twister (MT19937) for cross-platform
 *      reproducibility. Results will be identical on Linux, macOS, and other
 *      platforms. Note: Results will differ from original SecAct which uses
 *      platform-specific rand().
 *
 * =============================================================================
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

#include <string.h>
#include <math.h>
#include <stdlib.h>

#ifdef _OPENMP
    #include <omp.h>
#endif


/* =============================================================================
 * SECTION 1: RANDOM NUMBER GENERATION
 * =============================================================================
 *
 * GSL MT19937 is used for cross-platform reproducibility.
 * Native rand() implementation is kept for reference but unused.
 */

/**
 * Fisher-Yates shuffle using native rand() - matches original SecAct
 */
/*
static void shuffle_array_native(int *array, int n) {
    for (int i = 0; i < n - 1; i++) {
        int j = i + rand() / (RAND_MAX / (n - i) + 1);
        int tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}
*/

/**
 * Fisher-Yates shuffle using GSL RNG - cross-platform reproducible
 */
static void shuffle_array_gsl(gsl_rng *rng, int *array, int n) {
    for (int i = 0; i < n - 1; i++) {
        int j = i + (int)gsl_rng_uniform_int(rng, (unsigned long)(n - i));
        int tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}


/**
 * Generate permutation table using native rand() - for SecAct compatibility
 */

/*
static int* generate_permutation_table_native(int n, int nrand) {
    size_t table_size = (size_t)nrand * (size_t)n;
    int *table = (int*)malloc(table_size * sizeof(int));
    int *array = (int*)malloc((size_t)n * sizeof(int));

    // Same seed as original SecAct 
    srand(0);

    for (int i = 0; i < n; i++) {
        array[i] = i;
    }

    // Cumulative shuffle - matches original SecAct behavior
    for (int perm = 0; perm < nrand; perm++) {
        shuffle_array_native(array, n);
        memcpy(table + ((size_t)perm * n), array, n * sizeof(int));
    }

    free(array);
    return table;
}
*/

/**
 * Generate permutation table using GSL RNG - for cross-platform reproducibility
 */
static int* generate_permutation_table(int n, int nrand, unsigned long seed) {
    size_t table_size = (size_t)nrand * (size_t)n;
    int *table = (int*)malloc(table_size * sizeof(int));
    int *array = (int*)malloc((size_t)n * sizeof(int));

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    for (int i = 0; i < n; i++) {
        array[i] = i;
    }

    for (int perm = 0; perm < nrand; perm++) {
        shuffle_array_gsl(rng, array, n);
        memcpy(table + ((size_t)perm * n), array, n * sizeof(int));
    }

    gsl_rng_free(rng);
    free(array);
    return table;
}

/* Debug function to print GSL RNG values for Python comparison */
SEXP debug_gsl_rng(void) {
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 0);
    
    /* [1] Raw MT19937 values */
    Rprintf("[1] First 10 raw MT19937 values (seed=0):\n");
    for (int i = 0; i < 10; i++) {
        Rprintf("    %d: %lu\n", i, gsl_rng_get(rng));
    }
    
    /* [2] uniform_int(10) */
    gsl_rng_set(rng, 0);
    Rprintf("\n[2] uniform_int(10), first 10 values:\n    ");
    for (int i = 0; i < 10; i++) {
        Rprintf("%lu ", gsl_rng_uniform_int(rng, 10));
    }
    Rprintf("\n");
    
    /* [3] Shuffle [0..9] */
    gsl_rng_set(rng, 0);
    int arr[10] = {0,1,2,3,4,5,6,7,8,9};
    Rprintf("\n[3] Fisher-Yates shuffle of [0..9] (seed=0):\n");
    Rprintf("    Before: ");
    for (int i = 0; i < 10; i++) Rprintf("%d ", arr[i]);
    Rprintf("\n");
    
    for (int i = 0; i < 9; i++) {
        int j = i + (int)gsl_rng_uniform_int(rng, (unsigned long)(10 - i));
        int tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }
    Rprintf("    After:  ");
    for (int i = 0; i < 10; i++) Rprintf("%d ", arr[i]);
    Rprintf("\n");
    
    /* [4] Cumulative shuffles */
    gsl_rng_set(rng, 0);
    int arr2[10] = {0,1,2,3,4,5,6,7,8,9};
    Rprintf("\n[4] Cumulative shuffles of [0..9], 5 permutations (seed=0):\n");
    for (int perm = 0; perm < 5; perm++) {
        for (int i = 0; i < 9; i++) {
            int j = i + (int)gsl_rng_uniform_int(rng, (unsigned long)(10 - i));
            int tmp = arr2[j];
            arr2[j] = arr2[i];
            arr2[i] = tmp;
        }
        Rprintf("    Perm %d: ", perm);
        for (int i = 0; i < 10; i++) Rprintf("%d ", arr2[i]);
        Rprintf("\n");
    }
    
    gsl_rng_free(rng);
    return R_NilValue;
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
    mat->owner = 1;

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
/*
static void compute_projection_matrix(const gsl_matrix *Xt, double lambda,
                                       gsl_matrix *I_out, gsl_matrix *T_out) {
    // I = λI + X'X 
    gsl_matrix_set_identity(I_out);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_out);

    // I = (X'X + λI)^{-1} via Cholesky decomposition 
    gsl_linalg_cholesky_decomp(I_out);
    gsl_linalg_cholesky_invert(I_out);

    // T = (X'X + λI)^{-1} X' 
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_out, Xt, 0.0, T_out);
}
*/

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
/*
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
*/

/* =============================================================================
 * SECTION 4: LEGACY - .C INTERFACE (Y-PERMUTATION)
 * =============================================================================
 *
 * Original implementation using R's .C() interface.
 * Single-threaded with Y-row permutation strategy.
 * This version matches the original SecAct C code exactly.
 */

void ridgeReg(double *X_vec, double *Y_vec,
              int *n_ptr, int *p_ptr, int *m_ptr,
              double *lambda_ptr, double *nrand_ptr,
              double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec) {

    /* Extract dimensions - same as original */
    int n = *n_ptr;       /* Number of genes (rows) */
    int p = *p_ptr;       /* Number of proteins (columns of X) */
    int m = *m_ptr;       /* Number of samples (columns of Y) */
    int nrand = (int)*nrand_ptr;
    double lambda = *lambda_ptr;
    size_t total_elements = (size_t)p * (size_t)m;

    /* Wrap R vectors as GSL matrices - SAME DIMENSIONS AS ORIGINAL */
    gsl_matrix *X = wrap_r_vector_as_gsl_matrix(X_vec, (size_t)n, (size_t)p);
    gsl_matrix *Y = wrap_r_vector_as_gsl_matrix(Y_vec, (size_t)n, (size_t)m);
    gsl_matrix *beta = wrap_r_vector_as_gsl_matrix(beta_vec, (size_t)p, (size_t)m);
    gsl_matrix *aver_sq = wrap_r_vector_as_gsl_matrix(se_vec, (size_t)p, (size_t)m);
    gsl_matrix *zscore = wrap_r_vector_as_gsl_matrix(zscore_vec, (size_t)p, (size_t)m);
    gsl_matrix *pvalue = wrap_r_vector_as_gsl_matrix(pvalue_vec, (size_t)p, (size_t)m);

    /* Allocate working matrices */
    gsl_matrix *I = gsl_matrix_alloc((size_t)p, (size_t)p);
    gsl_matrix *T = gsl_matrix_alloc((size_t)p, (size_t)n);
    gsl_matrix *Y_rand = gsl_matrix_alloc((size_t)n, (size_t)m);
    gsl_matrix *beta_rand = gsl_matrix_alloc((size_t)p, (size_t)m);
    gsl_matrix *aver = gsl_matrix_alloc((size_t)p, (size_t)m);

    /* Compute beta = (X'X + lambda*I)^-1 * X' * Y */
    /* I = (X'X + lambda*I) */
    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);

    /* I = (X'X + lambda*I)^-1 */
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);

    /* T = (X'X + lambda*I)^-1 * X' */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);

    /* beta = (X'X + lambda*I)^-1 * X' * Y */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

    /* Permutation testing - use GSL RNG for cross-platform reproducibility */
    int *perm_table = generate_permutation_table(n, nrand, 0);
    gsl_matrix_set_zero(aver);
    gsl_matrix_set_zero(aver_sq);
    gsl_matrix_set_zero(pvalue);

    for (int i = 0; i < nrand; i++) {
        int *perm_idx = perm_table + ((size_t)i * (size_t)n);

        /* Create randomized Y by permuting rows */
        for (int j = 0; j < n; j++) {
            gsl_vector_const_view row = gsl_matrix_const_row(Y, perm_idx[j]);
            gsl_matrix_set_row(Y_rand, j, &row.vector);
        }

        /* Compute permuted coefficients */
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);

        /* P-value comparison */
        for (size_t j = 0; j < total_elements; j++) {
            if (fabs(beta_rand->data[j]) >= fabs(beta->data[j])) {
                pvalue->data[j]++;
            }
        }

        /* Accumulate for variance calculation */
        gsl_matrix_add(aver, beta_rand);
        gsl_matrix_mul_elements(beta_rand, beta_rand);
        gsl_matrix_add(aver_sq, beta_rand);
    }

    /* Finalize statistics - matching original exactly */
    gsl_matrix_scale(aver, 1.0 / nrand);
    gsl_matrix_scale(aver_sq, 1.0 / nrand);
    gsl_matrix_add_constant(pvalue, 1.0);
    int nrand1 = nrand + 1;
    gsl_matrix_scale(pvalue, 1.0 / nrand1);

    /* Compute z-score */
    gsl_matrix_memcpy(zscore, beta);
    gsl_matrix_sub(zscore, aver);
    gsl_matrix_mul_elements(aver, aver);
    gsl_matrix_sub(aver_sq, aver);
    for (size_t i = 0; i < total_elements; i++) {
        aver_sq->data[i] = sqrt(aver_sq->data[i]);
    }
    gsl_matrix_div_elements(zscore, aver_sq);

    /* Cleanup */
    free(perm_table);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Y_rand);
    gsl_matrix_free(beta_rand);
    gsl_matrix_free(aver);
    gsl_matrix_partial_free(X);
    gsl_matrix_partial_free(Y);
    gsl_matrix_partial_free(beta);
    gsl_matrix_partial_free(aver_sq);
    gsl_matrix_partial_free(zscore);
    gsl_matrix_partial_free(pvalue);
}


/* =============================================================================
 * SECTION 5: STYP - SINGLE-THREAD Y-PERMUTATION (.Call)
 * =============================================================================
 *
 * Uses same algorithm as legacy but with .Call interface.
 * R passes non-transposed matrices, C handles the layout conversion.
 */

SEXP ridgeReg_styp_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp) {

    /* Extract dimensions - R matrices are n×p and n×m */
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    int n = INTEGER(x_dim)[0];  /* Number of genes (rows) */
    int p = INTEGER(x_dim)[1];  /* Number of proteins (columns) */
    int m = INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];  /* Number of samples */
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

    memset(bv, 0, len * sizeof(double));
    memset(pv, 0, len * sizeof(double));
    memset(sv, 0, len * sizeof(double));
    memset(zv, 0, len * sizeof(double));

    /* 
     * R stores matrices in column-major order.
     * For an R matrix X (n×p), the data layout is: X[,1], X[,2], ..., X[,p]
     * 
     * If we create gsl_matrix_view_array(data, p, n), GSL interprets this as
     * a row-major p×n matrix, which gives us X^T.
     * 
     * So: Xt in GSL = X^T in mathematical terms
     */
    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);  /* X^T: p×n */
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);  /* Y^T: m×n */

    /* Allocate working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y = gsl_matrix_alloc(n, m);      /* Will hold Y (not transposed) */
    gsl_matrix *Y_rand = gsl_matrix_alloc(n, m);
    gsl_matrix *beta = gsl_matrix_alloc(p, m);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
    gsl_matrix *aver = gsl_matrix_alloc(p, m);
    gsl_matrix *aver_sq = gsl_matrix_alloc(p, m);

    /* Compute I = (X^T X + lambda*I)^{-1} */
    /* Xt is p×n, so Xt * Xt^T = X^T * X = p×p */
    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);

    /* T = I * Xt = (X^T X + lambda*I)^{-1} * X^T : p×n */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);

    /* Get Y in proper orientation: transpose Yt to get Y (n×m) */
    gsl_matrix_transpose_memcpy(Y, &Yt.matrix);

    /* beta = T * Y = p×n * n×m = p×m */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

    /* Copy beta to output */
    memcpy(bv, beta->data, len * sizeof(double));

    /* Permutation testing - use GSL RNG for cross-platform reproducibility */
    int *perm_table = generate_permutation_table(n, nrand, 0);
    gsl_matrix_set_zero(aver);
    gsl_matrix_set_zero(aver_sq);

    for (int i = 0; i < nrand; i++) {
        int *perm_idx = perm_table + ((size_t)i * n);

        /* Permute rows of Y */
        for (int j = 0; j < n; j++) {
            gsl_vector_const_view row = gsl_matrix_const_row(Y, perm_idx[j]);
            gsl_matrix_set_row(Y_rand, j, &row.vector);
        }

        /* Compute permuted beta */
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);

        /* Accumulate statistics */
        for (R_xlen_t k = 0; k < len; k++) {
            if (fabs(beta_rand->data[k]) >= fabs(beta->data[k])) pv[k]++;
        }
        gsl_matrix_add(aver, beta_rand);
        gsl_matrix_mul_elements(beta_rand, beta_rand);
        gsl_matrix_add(aver_sq, beta_rand);
    }

    /* Finalize statistics - matching original exactly */
    gsl_matrix_scale(aver, 1.0 / nrand);
    gsl_matrix_scale(aver_sq, 1.0 / nrand);

    /* Compute z-score and SE */
    for (R_xlen_t k = 0; k < len; k++) {
        pv[k] = (pv[k] + 1.0) / (nrand + 1.0);
        double mean = aver->data[k];
        double var = aver_sq->data[k] - mean * mean;
        sv[k] = sqrt(var);
        zv[k] = (sv[k] > 1e-12) ? (bv[k] - mean) / sv[k] : 0.0;
    }

    /* Cleanup */
    free(perm_table);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Y);
    gsl_matrix_free(Y_rand);
    gsl_matrix_free(beta);
    gsl_matrix_free(beta_rand);
    gsl_matrix_free(aver);
    gsl_matrix_free(aver_sq);

    SEXP result = create_result_list(beta_sexp, se_sexp, zscore_sexp, pvalue_sexp);
    UNPROTECT(4);
    return result;
}


/* =============================================================================
 * SECTION 6: STTP - SINGLE-THREAD T-PERMUTATION (.Call)
 * =============================================================================
 *
 * Uses T-column permutation (scatter approach) which can be more cache-friendly
 * for certain matrix sizes.
 */

SEXP ridgeReg_sttp_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp) {

    /* Extract dimensions */
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    int n = INTEGER(x_dim)[0];
    int p = INTEGER(x_dim)[1];
    int m = INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];
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

    memset(bv, 0, len * sizeof(double));
    memset(pv, 0, len * sizeof(double));
    memset(sv, 0, len * sizeof(double));
    memset(zv, 0, len * sizeof(double));

    /* Create GSL views - Xt = X^T (p×n), Yt = Y^T (m×n) */
    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);

    /* Allocate working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y = gsl_matrix_alloc(n, m);
    gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);   /* T transposed, original */
    gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);   /* T transposed, permuted */
    gsl_matrix *beta = gsl_matrix_alloc(p, m);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
    gsl_matrix *aver = gsl_matrix_alloc(p, m);
    gsl_matrix *aver_sq = gsl_matrix_alloc(p, m);

    /* Compute I = (X^T X + lambda*I)^{-1} */
    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);

    /* T = I * Xt : p×n */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);

    /* Get Y and compute beta */
    gsl_matrix_transpose_memcpy(Y, &Yt.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);
    memcpy(bv, beta->data, len * sizeof(double));

    /* Prepare T transpose for permutation */
    gsl_matrix_transpose_memcpy(Tt_orig, T);

    /* Permutation testing - use GSL RNG for cross-platform reproducibility */
    int *perm_table = generate_permutation_table(n, nrand, 0);
    gsl_matrix_set_zero(aver);
    gsl_matrix_set_zero(aver_sq);

    for (int i = 0; i < nrand; i++) {
        int *perm_idx = perm_table + ((size_t)i * n);

        /* Scatter permutation: Tt_perm[perm[j], :] = Tt_orig[j, :] */
        for (int j = 0; j < n; j++) {
            memcpy(Tt_perm->data + ((size_t)perm_idx[j] * p),
                   Tt_orig->data + ((size_t)j * p),
                   p * sizeof(double));
        }

        /* beta_rand = Tt_perm^T * Y = T_perm * Y */
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Tt_perm, Y, 0.0, beta_rand);

        /* Accumulate statistics */
        for (R_xlen_t k = 0; k < len; k++) {
            if (fabs(beta_rand->data[k]) >= fabs(beta->data[k])) pv[k]++;
        }
        gsl_matrix_add(aver, beta_rand);
        gsl_matrix_mul_elements(beta_rand, beta_rand);
        gsl_matrix_add(aver_sq, beta_rand);
    }

    /* Finalize statistics */
    gsl_matrix_scale(aver, 1.0 / nrand);
    gsl_matrix_scale(aver_sq, 1.0 / nrand);

    for (R_xlen_t k = 0; k < len; k++) {
        pv[k] = (pv[k] + 1.0) / (nrand + 1.0);
        double mean = aver->data[k];
        double var = aver_sq->data[k] - mean * mean;
        sv[k] = sqrt(var);
        zv[k] = (sv[k] > 1e-12) ? (bv[k] - mean) / sv[k] : 0.0;
    }

    /* Cleanup */
    free(perm_table);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Y);
    gsl_matrix_free(Tt_orig);
    gsl_matrix_free(Tt_perm);
    gsl_matrix_free(beta);
    gsl_matrix_free(beta_rand);
    gsl_matrix_free(aver);
    gsl_matrix_free(aver_sq);

    SEXP result = create_result_list(beta_sexp, se_sexp, zscore_sexp, pvalue_sexp);
    UNPROTECT(4);
    return result;
}


/* =============================================================================
 * SECTION 7: MTYP - MULTI-THREAD Y-PERMUTATION (OpenMP)
 * =============================================================================
 */

SEXP ridgeReg_mtyp_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp,
                            SEXP nrand_sexp, SEXP ncores_sexp) {

    /* Configure OpenMP threads */
    int ncores = asInteger(ncores_sexp);
    #ifdef _OPENMP
        if (ncores > 0) omp_set_num_threads(ncores);
    #endif

    /* Extract dimensions */
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    int n = INTEGER(x_dim)[0];
    int p = INTEGER(x_dim)[1];
    int m = INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];
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

    memset(bv, 0, len * sizeof(double));
    memset(pv, 0, len * sizeof(double));
    memset(sv, 0, len * sizeof(double));
    memset(zv, 0, len * sizeof(double));

    /* Create GSL views */
    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);

    /* Allocate shared working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y = gsl_matrix_alloc(n, m);
    gsl_matrix *beta = gsl_matrix_alloc(p, m);

    /* Compute I = (X^T X + lambda*I)^{-1} */
    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);

    /* T = I * Xt */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);

    /* Get Y and compute beta */
    gsl_matrix_transpose_memcpy(Y, &Yt.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);
    memcpy(bv, beta->data, len * sizeof(double));

    /* Prepare permutation table - use GSL RNG for cross-platform reproducibility */
    int *perm_table = generate_permutation_table(n, nrand, 0);

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
    double *thread_aver = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *thread_aver_sq = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *thread_counts = (double*)calloc((size_t)actual_threads * len, sizeof(double));

    /* Parallel permutation loop */
    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
            tid = omp_get_thread_num();
        #endif

        double *my_aver = &thread_aver[(size_t)tid * len];
        double *my_aver_sq = &thread_aver_sq[(size_t)tid * len];
        double *my_counts = &thread_counts[(size_t)tid * len];

        /* Thread-local working matrices */
        gsl_matrix *Y_rand = gsl_matrix_alloc(n, m);
        gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < nrand; i++) {
            int *perm_idx = perm_table + ((size_t)i * n);

            /* Permute rows of Y */
            for (int j = 0; j < n; j++) {
                gsl_vector_const_view row = gsl_matrix_const_row(Y, perm_idx[j]);
                gsl_matrix_set_row(Y_rand, j, &row.vector);
            }

            /* Compute permuted beta */
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);

            /* Accumulate statistics */
            for (R_xlen_t k = 0; k < len; k++) {
                if (fabs(beta_rand->data[k]) >= fabs(beta->data[k])) my_counts[k]++;
                my_aver[k] += beta_rand->data[k];
                my_aver_sq[k] += beta_rand->data[k] * beta_rand->data[k];
            }
        }

        gsl_matrix_free(Y_rand);
        gsl_matrix_free(beta_rand);
    }

    /* Reduce thread-local results */
    for (int t = 0; t < actual_threads; t++) {
        size_t offset = (size_t)t * len;
        for (R_xlen_t k = 0; k < len; k++) {
            sv[k] += thread_aver[offset + k];
            zv[k] += thread_aver_sq[offset + k];
            pv[k] += thread_counts[offset + k];
        }
    }

    /* Finalize statistics */
    for (R_xlen_t k = 0; k < len; k++) {
        double mean = sv[k] / nrand;
        double var = (zv[k] / nrand) - mean * mean;
        pv[k] = (pv[k] + 1.0) / (nrand + 1.0);
        zv[k] = (sqrt(var) > 1e-12) ? (bv[k] - mean) / sqrt(var) : 0.0;
        sv[k] = sqrt(var);
    }

    /* Cleanup */
    free(thread_aver);
    free(thread_aver_sq);
    free(thread_counts);
    free(perm_table);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Y);
    gsl_matrix_free(beta);

    SEXP result = create_result_list(beta_sexp, se_sexp, zscore_sexp, pvalue_sexp);
    UNPROTECT(4);
    return result;
}


/* =============================================================================
 * SECTION 8: MTTP - MULTI-THREAD T-PERMUTATION (OpenMP)
 * =============================================================================
 */

SEXP ridgeReg_mttp_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp,
                              SEXP nrand_sexp, SEXP ncores_sexp) {

    /* Configure OpenMP threads */
    int ncores = asInteger(ncores_sexp);
    #ifdef _OPENMP
        if (ncores > 0) omp_set_num_threads(ncores);
    #endif

    /* Extract dimensions */
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    int n = INTEGER(x_dim)[0];
    int p = INTEGER(x_dim)[1];
    int m = INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];
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

    memset(bv, 0, len * sizeof(double));
    memset(pv, 0, len * sizeof(double));
    memset(sv, 0, len * sizeof(double));
    memset(zv, 0, len * sizeof(double));

    /* Create GSL views */
    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);

    /* Allocate shared working matrices */
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y = gsl_matrix_alloc(n, m);
    gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
    gsl_matrix *beta = gsl_matrix_alloc(p, m);

    /* Compute I = (X^T X + lambda*I)^{-1} */
    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);

    /* T = I * Xt */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);

    /* Get Y and compute beta */
    gsl_matrix_transpose_memcpy(Y, &Yt.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);
    memcpy(bv, beta->data, len * sizeof(double));

    /* Prepare T transpose for permutation */
    gsl_matrix_transpose_memcpy(Tt_orig, T);

    /* Prepare permutation table - use GSL RNG for cross-platform reproducibility */
    int *perm_table = generate_permutation_table(n, nrand, 0);

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
    double *thread_aver = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *thread_aver_sq = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *thread_counts = (double*)calloc((size_t)actual_threads * len, sizeof(double));

    /* Parallel permutation loop */
    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
            tid = omp_get_thread_num();
        #endif

        double *my_aver = &thread_aver[(size_t)tid * len];
        double *my_aver_sq = &thread_aver_sq[(size_t)tid * len];
        double *my_counts = &thread_counts[(size_t)tid * len];

        /* Thread-local working matrices */
        gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);
        gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

        #pragma omp for schedule(dynamic, 16)
        for (int i = 0; i < nrand; i++) {
            int *perm_idx = perm_table + ((size_t)i * n);

            /* Scatter permutation */
            for (int j = 0; j < n; j++) {
                memcpy(Tt_perm->data + ((size_t)perm_idx[j] * p),
                       Tt_orig->data + ((size_t)j * p),
                       p * sizeof(double));
            }

            /* Compute permuted beta */
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Tt_perm, Y, 0.0, beta_rand);

            /* Accumulate statistics */
            for (R_xlen_t k = 0; k < len; k++) {
                if (fabs(beta_rand->data[k]) >= fabs(beta->data[k])) my_counts[k]++;
                my_aver[k] += beta_rand->data[k];
                my_aver_sq[k] += beta_rand->data[k] * beta_rand->data[k];
            }
        }

        gsl_matrix_free(Tt_perm);
        gsl_matrix_free(beta_rand);
    }

    /* Reduce thread-local results */
    for (int t = 0; t < actual_threads; t++) {
        size_t offset = (size_t)t * len;
        for (R_xlen_t k = 0; k < len; k++) {
            sv[k] += thread_aver[offset + k];
            zv[k] += thread_aver_sq[offset + k];
            pv[k] += thread_counts[offset + k];
        }
    }

    /* Finalize statistics */
    for (R_xlen_t k = 0; k < len; k++) {
        double mean = sv[k] / nrand;
        double var = (zv[k] / nrand) - mean * mean;
        pv[k] = (pv[k] + 1.0) / (nrand + 1.0);
        zv[k] = (sqrt(var) > 1e-12) ? (bv[k] - mean) / sqrt(var) : 0.0;
        sv[k] = sqrt(var);
    }

    /* Cleanup */
    free(thread_aver);
    free(thread_aver_sq);
    free(thread_counts);
    free(perm_table);
    gsl_matrix_free(Tt_orig);
    gsl_matrix_free(I);
    gsl_matrix_free(T);
    gsl_matrix_free(Y);
    gsl_matrix_free(beta);

    SEXP result = create_result_list(beta_sexp, se_sexp, zscore_sexp, pvalue_sexp);
    UNPROTECT(4);
    return result;
}


/* =============================================================================
 * SECTION 9: R REGISTRATION
 * =============================================================================
 */

/* Add the debug function declaration at the top of the file or in a header */
SEXP debug_gsl_rng(void);

static const R_CMethodDef cMethods[] = {
    {"ridgeReg", (DL_FUNC) &ridgeReg, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef callMethods[] = {
    {"ridgeReg_styp_interface",      (DL_FUNC) &ridgeReg_styp_interface,      4},
    {"ridgeReg_sttp_interface",      (DL_FUNC) &ridgeReg_sttp_interface,      4},
    {"ridgeReg_mtyp_interface",      (DL_FUNC) &ridgeReg_mtyp_interface,      5},
    {"ridgeReg_mttp_interface",      (DL_FUNC) &ridgeReg_mttp_interface,      5},
    {"debug_gsl_rng",                (DL_FUNC) &debug_gsl_rng,                0},
    {NULL, NULL, 0}
};

void R_init_RidgeR(DllInfo *dll) {
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
