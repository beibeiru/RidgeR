#include <R.h>
#include <Rinternals.h> // Required for .Call
#include <R_ext/Rdynload.h> // For registration
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>  // For GSL RNG (cross-platform reproducibility)
#include <string.h> // for memcpy
#include <time.h>   // for seed
#include <math.h>   // for sqrt, fabs

#ifdef _OPENMP
  #undef match        // R's Rinternals.h defines match as Rf_match, which
                      // conflicts with OpenMP's "match" clause in omp.h
  #include <omp.h>
  #define match Rf_match
#endif

// --- TUNING ---
// Number of samples (columns) a thread processes at one time.
// 64 is small enough to fit in cache, large enough for BLAS efficiency.
// This keeps memory usage CONSTANT per thread, regardless of m.
#define SAMPLE_STRIP_SIZE 64

// Batch of permutations.
#define PERM_BATCH_SIZE 64

// GSL MT19937 max value (2^32 - 1)
#define GSL_MT19937_MAX 0xFFFFFFFFUL

/* =========================================================
   RNG METHOD CONSTANTS
   0 = srand/rand (C stdlib, platform-dependent, matches old behavior)
   1 = GSL MT19937 (cross-platform reproducible, matches SecActPy GSLRNG)
   ========================================================= */
#define RNG_SRAND 0
#define RNG_GSL   1

/* =========================================================
   HELPER FUNCTIONS
   ========================================================= */

gsl_matrix *RVectorObject_to_gsl_matrix(double *vec, size_t nr, size_t nc)
{
  gsl_block *b = (gsl_block*)malloc(sizeof(gsl_block));
  gsl_matrix *r = (gsl_matrix*)malloc(sizeof(gsl_matrix));
  r->size1 = nr;
  r->tda = r->size2 = nc;
  r->owner = 1;
  b->data = r->data = vec;
  r->block = b;
  b->size = r->size1 * r->size2;
  return r;
}

void gsl_matrix_partial_free(gsl_matrix *x)
{
  if(x) {
    if(x->block) free(x->block);
    free(x);
  }
}

/* Fisher-Yates shuffle using C stdlib rand().
   Uses the formula: j = i + rand() / (RAND_MAX / (n - i) + 1)
   Platform-dependent (glibc rand() differs from macOS). */
void shuffle_srand(int array[], const int n)
{
  int i, j;
  double t;
  for (i = 0; i < n-1; i++) {
    j = i + rand() / (RAND_MAX / (n - i) + 1);
    t = array[j];
    array[j] = array[i];
    array[i] = t;
  }
}

/* Fisher-Yates shuffle using GSL MT19937.
   Same formula as shuffle_srand but with gsl_rng_get() and GSL_MT19937_MAX.
   Cross-platform reproducible. Matches SecActPy's GSLRNG.shuffle_inplace(). */
void shuffle_gsl(gsl_rng *rng, int array[], const int n)
{
  int i, j;
  double t;
  for (i = 0; i < n-1; i++) {
    unsigned long r = gsl_rng_get(rng);
    j = i + (int)(r / (GSL_MT19937_MAX / (unsigned long)(n - i) + 1));
    t = array[j];
    array[j] = array[i];
    array[i] = t;
  }
}

/* =========================================================
   OLD VERSION (updated with rng_method parameter)
   ========================================================= */
void ridgeReg(
  double *X_vec, double *Y_vec,
  int *n_pt, int *p_pt, int *m_pt,
  double *lambda_pt, double *nrand_pt,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int *rng_method_pt
)
{
  gsl_matrix *X, *Y, *I, *T, *beta, *Y_rand, *beta_rand, *aver, *aver_sq, *zscore, *pvalue;
  int n = *n_pt, p = *p_pt, m = *m_pt, nrand = (int)*nrand_pt;
  double lambda = *lambda_pt;
  int rng_method = *rng_method_pt;
  int *array_index, i, j;

  X = RVectorObject_to_gsl_matrix(X_vec, n, p);
  Y = RVectorObject_to_gsl_matrix(Y_vec, n, m);
  beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);
  aver_sq = RVectorObject_to_gsl_matrix(se_vec, p, m);
  zscore = RVectorObject_to_gsl_matrix(zscore_vec, p, m);
  pvalue = RVectorObject_to_gsl_matrix(pvalue_vec, p, m);

  I = gsl_matrix_alloc(p, p);
  T = gsl_matrix_alloc(p, n);
  Y_rand = gsl_matrix_alloc(n, m);
  beta_rand = gsl_matrix_alloc(p, m);
  aver = gsl_matrix_alloc(p, m);
  array_index = (int*)malloc(n*sizeof(int));

  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

  for(i=0;i<n;i++) array_index[i] = i;

  gsl_matrix_set_zero(aver);
  gsl_matrix_set_zero(aver_sq);
  gsl_matrix_set_zero(pvalue);

  /* Initialize RNG based on method */
  gsl_rng *gsl_rng_obj = NULL;
  if (rng_method == RNG_GSL) {
    gsl_rng_obj = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gsl_rng_obj, 0);  /* seed=0 -> GSL uses 4357 internally */
  } else {
    srand(0);
  }

  for(i=0;i<nrand;i++)
  {
    if (rng_method == RNG_GSL) {
      shuffle_gsl(gsl_rng_obj, array_index, n);
    } else {
      shuffle_srand(array_index, n);
    }

    for(j=0;j<n;j++){
      gsl_vector_const_view t = gsl_matrix_const_row(Y, array_index[j]);
      gsl_matrix_set_row(Y_rand, j, &t.vector);
    }

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);

    for(j=0; j< pvalue->size1 * pvalue->size2; j++) {
      if(fabs(beta_rand->data[j]) >= fabs(beta->data[j])) pvalue->data[j]++;
    }

    gsl_matrix_add(aver, beta_rand);
    gsl_matrix_mul_elements(beta_rand, beta_rand);
    gsl_matrix_add(aver_sq, beta_rand);
  }

  gsl_matrix_scale(aver, 1.0/ nrand);
  gsl_matrix_scale(aver_sq, 1.0/ nrand);

  gsl_matrix_add_constant(pvalue, 1.0);
  gsl_matrix_scale(pvalue, 1.0/ (nrand+1));

  gsl_matrix_memcpy(zscore, beta);
  gsl_matrix_sub(zscore, aver);

  gsl_matrix_mul_elements(aver, aver);
  gsl_matrix_sub(aver_sq, aver);

  for(i=0;i< aver_sq->size1 * aver_sq->size2; i++) {
    aver_sq->data[i] = sqrt(aver_sq->data[i]);
  }

  gsl_matrix_div_elements(zscore, aver_sq);

  if (gsl_rng_obj) gsl_rng_free(gsl_rng_obj);
  gsl_matrix_free(I);
  gsl_matrix_free(T);
  gsl_matrix_free(Y_rand);
  gsl_matrix_free(beta_rand);
  gsl_matrix_free(aver);
  free(array_index);
  gsl_matrix_partial_free(beta);
  gsl_matrix_partial_free(aver_sq);
  gsl_matrix_partial_free(zscore);
  gsl_matrix_partial_free(pvalue);
}

/* =========================================================
   NEW OPTIMIZED FUNCTION (Long Vector Support)
   ========================================================= */

// Internal logic for the new function
void ridgeRegFast_core(
  double *X_ptr, double *Y_ptr, // Raw R pointers
  size_t n, size_t p, size_t m, // n=genes, p=proteins, m=samples
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int num_threads, int rng_method
)
{
  if (num_threads > 0) {
       #ifdef _OPENMP
           omp_set_num_threads(num_threads);
       #endif
  }

  // --- THE TRICK ---
  // R stores X as (n x p) Column-Major.
  // GSL reads Row-Major.
  // If we map X_ptr as (p x n), GSL sees exactly t(X).
  // This is perfect because the math usually requires X' anyway.

  // We want X_gsl to be (p x n) and Y_gsl to be (m x n) for efficient row slicing
  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n); // Effectively t(X)
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n); // Effectively t(Y)
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  // Allocations
  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);

  // 1. Projection: T = (X'X + lambda*I)^-1 * X'
  // Note: Since Xt is already transposed view, we use CblasNoTrans where we used Trans before
  // Math: X_real = X (n x p). Xt = X' (p x n).
  // We need (Xt * X)^-1.
  // In GSL with our view: Xt * Xt'

  gsl_matrix_set_identity(I_mat);
  // Xt (p x n) * Xt^T (n x p) -> (p x p)
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
  gsl_linalg_cholesky_decomp(I_mat);
  gsl_linalg_cholesky_invert(I_mat);

  // T = Inv * Xt
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);

  // 2. Real Beta = T * Y_real
  // Y_real is (n x m). Yt view is (m x n).
  // Beta = T (p x n) * Y_real (n x m).
  // Using Yt view: T * Yt^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  // 3. Permutations - generate table using selected RNG
  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  for(size_t k=0; k<n; k++) temp_idx[k] = (int)k;

  gsl_rng *gsl_rng_obj = NULL;
  if (rng_method == RNG_GSL) {
    gsl_rng_obj = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gsl_rng_obj, 0);  /* seed=0 -> GSL uses 4357 internally */
    for(int i_rand = 0; i_rand < nrand; i_rand++) {
      shuffle_gsl(gsl_rng_obj, temp_idx, (int)n);
      memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
    }
    gsl_rng_free(gsl_rng_obj);
  } else {
    srand(0);
    for(int i_rand = 0; i_rand < nrand; i_rand++) {
      shuffle_srand(temp_idx, (int)n);
      memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
    }
  }
  free(temp_idx);

  // Init outputs
  R_xlen_t total_elements = (R_xlen_t)p * (R_xlen_t)m;
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = 0.0; se_vec[i] = 0.0; zscore_vec[i] = 0.0;
  }

  // 4. Parallel Loop
  #pragma omp parallel
  {
    // Buffers
    gsl_matrix *Y_block = gsl_matrix_alloc(n, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    gsl_matrix *Beta_block = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);

    // Loop over SAMPLES (which are Rows in our Yt view!)
    // This allows efficient row slicing which GSL loves.
    #pragma omp for schedule(dynamic)
    for(size_t samp_start = 0; samp_start < m; samp_start += SAMPLE_STRIP_SIZE)
    {
      size_t current_samples = (samp_start + SAMPLE_STRIP_SIZE > m) ? (m - samp_start) : SAMPLE_STRIP_SIZE;

      for(int b_start = 0; b_start < nrand; b_start += PERM_BATCH_SIZE)
      {
        int current_perms = (b_start + PERM_BATCH_SIZE > nrand) ? (nrand - b_start) : PERM_BATCH_SIZE;

        // Construct Y_block
        // Yt is (m samples x n genes). We take specific rows (samples) and permute columns (genes).
        for(int i_perm = 0; i_perm < current_perms; i_perm++) {
           int *p_idx = &perm_table[(size_t)(b_start + i_perm) * n];

           for(size_t s_local = 0; s_local < current_samples; s_local++) {
             // We need to construct columns in Y_block from rows in Yt
             // This is a transposing copy, but block size is small (cache friendly).

             // Yt row index
             size_t global_samp_idx = samp_start + s_local;
             double *src_row = gsl_matrix_ptr(Yt, global_samp_idx, 0);

             // Destination column in Y_block
             size_t dest_col = (i_perm * current_samples) + s_local;

             for(size_t gene = 0; gene < n; gene++) {
               // Apply permutation here.
               // Standard logic: permute rows of Y.
               // Here Yt is (Samples x Genes). We want to permute Genes.
               // So we pick src_row[ p_idx[gene] ]
               double val = src_row[ p_idx[gene] ];
               gsl_matrix_set(Y_block, gene, dest_col, val);
             }
           }
        }

        // Multiply: Beta = T * Y_block
        gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block, 0, 0, n, current_perms * current_samples);
        gsl_matrix_view B_sub = gsl_matrix_submatrix(Beta_block, 0, 0, p, current_perms * current_samples);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);

        // Stats
        for(int i_perm = 0; i_perm < current_perms; i_perm++) {
          for(size_t r = 0; r < p; r++) {
            for(size_t s_local = 0; s_local < current_samples; s_local++) {

              R_xlen_t global_idx = (R_xlen_t)r * m + (samp_start + s_local);

              double b_rnd = gsl_matrix_get(&B_sub.matrix, r, ((size_t)i_perm * current_samples) + s_local);
              double b_obs = beta_vec[global_idx];

              if(fabs(b_rnd) >= fabs(b_obs)) pvalue_vec[global_idx] += 1.0;
              zscore_vec[global_idx] += b_rnd;
              se_vec[global_idx]     += (b_rnd * b_rnd);
            }
          }
        }
      }
    }
    gsl_matrix_free(Y_block); gsl_matrix_free(Beta_block);
  }

  free(perm_table);
  // Finalize (Same as before)
  double inv_nrand = 1.0 / nrand;
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean_rand = zscore_vec[i] * inv_nrand;
    double var_rand = (se_vec[i] * inv_nrand) - (mean_rand * mean_rand);
    double se_rand = sqrt(var_rand > 0 ? var_rand : 0);
    se_vec[i] = se_rand;
    zscore_vec[i] = (se_rand > 1e-12) ? (beta_vec[i] - mean_rand) / se_rand : 0.0;
  }

  gsl_matrix_free(I_mat); gsl_matrix_free(T);
  gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta);
}


// .Call Interface (Required for Long Vectors)
SEXP ridgeRegFast_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp, SEXP ncores_sexp, SEXP rng_method_sexp)
{
    // INPUTS ARE NOW UN-TRANSPOSED
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0]; // Genes
    size_t p = (size_t)INTEGER(x_dim)[1]; // Proteins

    SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
    // n from Y must match n from X
    size_t m = (size_t)INTEGER(y_dim)[1]; // Samples

    R_xlen_t total_len = (R_xlen_t)p * (R_xlen_t)m;

    SEXP beta_sexp   = PROTECT(allocVector(REALSXP, total_len));
    SEXP se_sexp     = PROTECT(allocVector(REALSXP, total_len));
    SEXP zscore_sexp = PROTECT(allocVector(REALSXP, total_len));
    SEXP pvalue_sexp = PROTECT(allocVector(REALSXP, total_len));

    ridgeRegFast_core(
        REAL(X_sexp), REAL(Y_sexp), n, p, m,
        asReal(lambda_sexp), asInteger(nrand_sexp),
        REAL(beta_sexp), REAL(se_sexp), REAL(zscore_sexp), REAL(pvalue_sexp),
        asInteger(ncores_sexp), asInteger(rng_method_sexp)
    );

    SEXP res_list = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(res_list, 0, beta_sexp);
    SET_VECTOR_ELT(res_list, 1, se_sexp);
    SET_VECTOR_ELT(res_list, 2, zscore_sexp);
    SET_VECTOR_ELT(res_list, 3, pvalue_sexp);

    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("beta"));
    SET_STRING_ELT(names, 1, mkChar("se"));
    SET_STRING_ELT(names, 2, mkChar("zscore"));
    SET_STRING_ELT(names, 3, mkChar("pvalue"));
    setAttrib(res_list, R_NamesSymbol, names);

    UNPROTECT(6);
    return res_list;
}

/* =========================================================
   GENERATE PERMUTATION TABLE (.Call)
   Returns integer matrix (nrand x n) of forward permutation indices.
   Uses cumulative Fisher-Yates (same RNG sequence as ridgeReg/ridgeRegFast).
   Enables the naive R variant to use identical RNG sequences.
   ========================================================= */
SEXP generate_perm_table(SEXP n_sexp, SEXP nrand_sexp, SEXP rng_method_sexp)
{
    int n = asInteger(n_sexp);
    int nrand = asInteger(nrand_sexp);
    int rng_method = asInteger(rng_method_sexp);

    /* Allocate R integer matrix (nrand x n) */
    SEXP result = PROTECT(allocMatrix(INTSXP, nrand, n));
    int *result_ptr = INTEGER(result);

    /* Working array for cumulative shuffle */
    int *temp_idx = (int*)malloc(n * sizeof(int));
    for (int k = 0; k < n; k++) temp_idx[k] = k;

    /* Initialize RNG */
    gsl_rng *gsl_rng_obj = NULL;
    if (rng_method == RNG_GSL) {
        gsl_rng_obj = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gsl_rng_obj, 0);
    } else {
        srand(0);
    }

    /* Generate permutations (cumulative Fisher-Yates) */
    for (int i = 0; i < nrand; i++) {
        if (rng_method == RNG_GSL) {
            shuffle_gsl(gsl_rng_obj, temp_idx, n);
        } else {
            shuffle_srand(temp_idx, n);
        }
        /* Store row i of the matrix (column-major: element [i, j] at j*nrand + i) */
        for (int j = 0; j < n; j++) {
            result_ptr[(size_t)j * nrand + i] = temp_idx[j];
        }
    }

    if (gsl_rng_obj) gsl_rng_free(gsl_rng_obj);
    free(temp_idx);

    UNPROTECT(1);
    return result;
}


/* =========================================================
   T-COLUMN PERMUTATION VARIANT (serves both st and mt)
   ========================================================= */

void ridgeRegFastTcol_core(
  double *X_ptr, double *Y_ptr,
  size_t n, size_t p, size_t m,
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int num_threads, int rng_method
)
{
  if (num_threads > 0) {
       #ifdef _OPENMP
           omp_set_num_threads(num_threads);
       #endif
  }

  /* Same transposed-view trick as ridgeRegFast_core */
  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n);
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);

  /* 1. Projection: T = (X'X + lambda*I)^-1 * X' */
  gsl_matrix_set_identity(I_mat);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
  gsl_linalg_cholesky_decomp(I_mat);
  gsl_linalg_cholesky_invert(I_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);

  /* 2. Observed beta = T * Y' (using Yt transposed view) */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  /* 3. Generate INVERSE permutation table for T-column permutation
     Forward shuffle temp_idx, then compute inv_perm[temp_idx[j]] = j */
  int *inv_perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  for (size_t k = 0; k < n; k++) temp_idx[k] = (int)k;

  gsl_rng *gsl_rng_obj = NULL;
  if (rng_method == RNG_GSL) {
    gsl_rng_obj = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gsl_rng_obj, 0);
    for (int i_rand = 0; i_rand < nrand; i_rand++) {
      shuffle_gsl(gsl_rng_obj, temp_idx, (int)n);
      /* Compute inverse permutation */
      int *inv_row = &inv_perm_table[(size_t)i_rand * n];
      for (size_t j = 0; j < n; j++) {
        inv_row[temp_idx[j]] = (int)j;
      }
    }
    gsl_rng_free(gsl_rng_obj);
  } else {
    srand(0);
    for (int i_rand = 0; i_rand < nrand; i_rand++) {
      shuffle_srand(temp_idx, (int)n);
      int *inv_row = &inv_perm_table[(size_t)i_rand * n];
      for (size_t j = 0; j < n; j++) {
        inv_row[temp_idx[j]] = (int)j;
      }
    }
  }
  free(temp_idx);

  /* Init outputs */
  R_xlen_t total_elements = (R_xlen_t)p * (R_xlen_t)m;
  #pragma omp parallel for schedule(static)
  for (R_xlen_t i = 0; i < total_elements; i++) {
    pvalue_vec[i] = 0.0; se_vec[i] = 0.0; zscore_vec[i] = 0.0;
  }

  /* 4. Parallel loop over permutations (T-column approach)
     Each thread: construct T_perm (p x n) from inv_perm, then process sample strips */
  #pragma omp parallel
  {
    /* Per-thread buffers */
    gsl_matrix *T_perm = gsl_matrix_alloc(p, n);
    gsl_matrix *beta_strip = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE);

    #pragma omp for schedule(dynamic)
    for (int i_rand = 0; i_rand < nrand; i_rand++)
    {
      int *inv_perm = &inv_perm_table[(size_t)i_rand * n];

      /* Construct T_perm: T_perm[:, j] = T[:, inv_perm[j]]
         i.e. permute columns of T */
      for (size_t col = 0; col < n; col++) {
        size_t src_col = (size_t)inv_perm[col];
        /* Copy column src_col of T to column col of T_perm */
        for (size_t row = 0; row < p; row++) {
          gsl_matrix_set(T_perm, row, col, gsl_matrix_get(T, row, src_col));
        }
      }

      /* Process sample strips */
      for (size_t samp_start = 0; samp_start < m; samp_start += SAMPLE_STRIP_SIZE)
      {
        size_t current_samples = (samp_start + SAMPLE_STRIP_SIZE > m)
                                  ? (m - samp_start) : SAMPLE_STRIP_SIZE;

        /* beta_strip = T_perm * Y_strip
           Y is stored as Yt (m x n). Y_strip is a sub-view of rows samp_start..+current_samples.
           We need T_perm (p x n) * Y_strip^T where Y_strip is (current_samples x n).
           Using Yt submatrix: */
        gsl_matrix_view Y_sub = gsl_matrix_submatrix(Yt, samp_start, 0, current_samples, n);
        gsl_matrix_view B_sub = gsl_matrix_submatrix(beta_strip, 0, 0, p, current_samples);

        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T_perm, &Y_sub.matrix, 0.0, &B_sub.matrix);

        /* Accumulate statistics (with atomics since multiple permutations
           may write to the same output elements) */
        for (size_t r = 0; r < p; r++) {
          for (size_t s_local = 0; s_local < current_samples; s_local++) {
            R_xlen_t global_idx = (R_xlen_t)r * m + (samp_start + s_local);
            double b_rnd = gsl_matrix_get(&B_sub.matrix, r, s_local);
            double b_obs = beta_vec[global_idx];

            #pragma omp atomic
            pvalue_vec[global_idx] += (fabs(b_rnd) >= fabs(b_obs)) ? 1.0 : 0.0;
            #pragma omp atomic
            zscore_vec[global_idx] += b_rnd;
            #pragma omp atomic
            se_vec[global_idx] += (b_rnd * b_rnd);
          }
        }
      }
    }
    gsl_matrix_free(T_perm);
    gsl_matrix_free(beta_strip);
  }

  free(inv_perm_table);

  /* 5. Finalize (same formula as ridgeRegFast_core) */
  double inv_nrand = 1.0 / nrand;
  #pragma omp parallel for schedule(static)
  for (R_xlen_t i = 0; i < total_elements; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean_rand = zscore_vec[i] * inv_nrand;
    double var_rand = (se_vec[i] * inv_nrand) - (mean_rand * mean_rand);
    double se_rand = sqrt(var_rand > 0 ? var_rand : 0);
    se_vec[i] = se_rand;
    zscore_vec[i] = (se_rand > 1e-12) ? (beta_vec[i] - mean_rand) / se_rand : 0.0;
  }

  gsl_matrix_free(I_mat); gsl_matrix_free(T);
  gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta);
}


// .Call Interface for T-column variant
SEXP ridgeRegFastTcol_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp, SEXP ncores_sexp, SEXP rng_method_sexp)
{
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0];
    size_t p = (size_t)INTEGER(x_dim)[1];

    SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
    size_t m = (size_t)INTEGER(y_dim)[1];

    R_xlen_t total_len = (R_xlen_t)p * (R_xlen_t)m;

    SEXP beta_sexp   = PROTECT(allocVector(REALSXP, total_len));
    SEXP se_sexp     = PROTECT(allocVector(REALSXP, total_len));
    SEXP zscore_sexp = PROTECT(allocVector(REALSXP, total_len));
    SEXP pvalue_sexp = PROTECT(allocVector(REALSXP, total_len));

    ridgeRegFastTcol_core(
        REAL(X_sexp), REAL(Y_sexp), n, p, m,
        asReal(lambda_sexp), asInteger(nrand_sexp),
        REAL(beta_sexp), REAL(se_sexp), REAL(zscore_sexp), REAL(pvalue_sexp),
        asInteger(ncores_sexp), asInteger(rng_method_sexp)
    );

    SEXP res_list = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(res_list, 0, beta_sexp);
    SET_VECTOR_ELT(res_list, 1, se_sexp);
    SET_VECTOR_ELT(res_list, 2, zscore_sexp);
    SET_VECTOR_ELT(res_list, 3, pvalue_sexp);

    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("beta"));
    SET_STRING_ELT(names, 1, mkChar("se"));
    SET_STRING_ELT(names, 2, mkChar("zscore"));
    SET_STRING_ELT(names, 3, mkChar("pvalue"));
    setAttrib(res_list, R_NamesSymbol, names);

    UNPROTECT(6);
    return res_list;
}


// Registration
static const R_CMethodDef cMethods[] = {
    {"ridgeReg", (DL_FUNC) &ridgeReg, 12},  /* 11 original + rng_method */
    {NULL, NULL, 0}
};

static const R_CallMethodDef callMethods[] = {
    {"ridgeRegFast_interface",     (DL_FUNC) &ridgeRegFast_interface,     6},
    {"ridgeRegFastTcol_interface", (DL_FUNC) &ridgeRegFastTcol_interface, 6},
    {"generate_perm_table",        (DL_FUNC) &generate_perm_table,        3},
    {NULL, NULL, 0}
};

void R_init_RidgeR(DllInfo *dll) {
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
