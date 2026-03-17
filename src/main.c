#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_version.h>

/* gsl_linalg_cholesky_decomp was deprecated in GSL 2.6;
   use gsl_linalg_cholesky_decomp1 when available.           */
#if GSL_MAJOR_VERSION > 2 || (GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 6)
#define CHOLESKY_DECOMP gsl_linalg_cholesky_decomp1
#else
#define CHOLESKY_DECOMP gsl_linalg_cholesky_decomp
#endif
#include <string.h>
#include <time.h>
#include <math.h>

/* =========================================================
   OPENMP — disabled on macOS+vecLib at compile time via
   configure. On vecLib, GCD already threads BLAS internally;
   adding OpenMP causes oversubscription and system slowdown.
   ========================================================= */
#ifdef _OPENMP
  #undef match      /* Rinternals.h defines match as Rf_match, conflicts with omp.h */
  #include <omp.h>
  #define match Rf_match
#endif

/* --- TUNING (multi-threaded paths only) ---
   SAMPLE_STRIP_SIZE: columns per thread per iteration.
   PERM_BATCH_SIZE  : permutations batched per iteration (Yrow.mt).
   These only apply to the parallel code paths. Single-threaded
   paths bypass batching entirely and use one dgemm per permutation. */
#define SAMPLE_STRIP_SIZE 64
#define PERM_BATCH_SIZE   64

/* GSL MT19937 max (2^32 - 1) */
#define GSL_MT19937_MAX 0xFFFFFFFFUL

/* RNG method constants */
#define RNG_SRAND 0   /* C stdlib rand()  — platform-dependent       */
#define RNG_GSL   1   /* GSL MT19937      — cross-platform, SecActPy  */


/* =========================================================
   HELPER FUNCTIONS
   ========================================================= */

/* Wrap a pre-allocated double* as a GSL matrix view (zero-copy).
   Use gsl_matrix_partial_free() to release the wrapper only. */
gsl_matrix *RVectorObject_to_gsl_matrix(double *vec, size_t nr, size_t nc)
{
  gsl_block  *b = (gsl_block*)malloc(sizeof(gsl_block));
  gsl_matrix *r = (gsl_matrix*)malloc(sizeof(gsl_matrix));
  r->size1 = nr;
  r->tda   = r->size2 = nc;
  r->owner = 1;
  b->data  = r->data = vec;
  r->block = b;
  b->size  = nr * nc;
  return r;
}

/* Free only the GSL wrapper structs, not the underlying data buffer. */
void gsl_matrix_partial_free(gsl_matrix *x)
{
  if (x) { if (x->block) free(x->block); free(x); }
}

/* Fisher-Yates shuffle — C stdlib rand() (platform-dependent). */
void shuffle_srand(int array[], const int n)
{
  int i, j, t;
  for (i = 0; i < n-1; i++) {
    j = i + rand() / (RAND_MAX / (n - i) + 1);
    t = array[j]; array[j] = array[i]; array[i] = t;
  }
}

/* Fisher-Yates shuffle — GSL MT19937 (cross-platform, matches SecActPy). */
void shuffle_gsl(gsl_rng *rng, int array[], const int n)
{
  int i, j, t;
  for (i = 0; i < n-1; i++) {
    unsigned long r = gsl_rng_get(rng);
    j = i + (int)(r / (GSL_MT19937_MAX / (unsigned long)(n - i) + 1));
    t = array[j]; array[j] = array[i]; array[i] = t;
  }
}

/* Build the forward permutation table (nrand x n).
   Shared by both Yrow and Tcol cores to guarantee identical RNG streams. */
static int *build_perm_table(int nrand, size_t n, int rng_method)
{
  int *table    = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  for (size_t k = 0; k < n; k++) temp_idx[k] = (int)k;

  gsl_rng *gsl_rng_obj = NULL;
  if (rng_method == RNG_GSL) {
    gsl_rng_obj = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gsl_rng_obj, 0);
    for (int i = 0; i < nrand; i++) {
      shuffle_gsl(gsl_rng_obj, temp_idx, (int)n);
      memcpy(&table[(size_t)i * n], temp_idx, n * sizeof(int));
    }
    gsl_rng_free(gsl_rng_obj);
  } else {
    srand(0);
    for (int i = 0; i < nrand; i++) {
      shuffle_srand(temp_idx, (int)n);
      memcpy(&table[(size_t)i * n], temp_idx, n * sizeof(int));
    }
  }
  free(temp_idx);
  return table;
}


/* Finalize permutation statistics: convert accumulated sums into
   p-value, SE, and z-score. Shared by both Yrow and Tcol cores. */
static void finalize_stats(
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  R_xlen_t total, int nrand)
{
  double inv_nrand = 1.0 / nrand;
  for (R_xlen_t i = 0; i < total; i++) {
    pvalue_vec[i]  = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean_r  = zscore_vec[i] * inv_nrand;
    double var_r   = (se_vec[i] * inv_nrand) - mean_r * mean_r;
    double se_r    = sqrt(var_r > 0 ? var_r : 0);
    se_vec[i]      = se_r;
    zscore_vec[i]  = (se_r > 1e-12) ? (beta_vec[i] - mean_r) / se_r : 0.0;
  }
}

/* Accumulate permutation statistics (serial, no atomics).
   Used by both Yrow and Tcol single-threaded paths. */
static inline void accumulate_stats_serial(
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  gsl_matrix *beta_rand, size_t p, size_t m)
{
  for (size_t r = 0; r < p; r++) {
    for (size_t s = 0; s < m; s++) {
      R_xlen_t idx   = (R_xlen_t)r * m + s;
      double   b_rnd = beta_rand->data[r * beta_rand->tda + s];
      double   b_obs = beta_vec[idx];
      if (fabs(b_rnd) >= fabs(b_obs)) pvalue_vec[idx] += 1.0;
      zscore_vec[idx] += b_rnd;
      se_vec[idx]     += b_rnd * b_rnd;
    }
  }
}

/* Compute projection matrix: T = (X'X + lambda*I)^-1 * X'.
   Xt is the p-by-n transposed design matrix (GSL row-major).
   I_mat (p x p) and T (p x n) must be pre-allocated. */
static void compute_projection(gsl_matrix *Xt, gsl_matrix *I_mat, gsl_matrix *T,
                                size_t p, double lambda)
{
  gsl_matrix_set_identity(I_mat);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
  CHOLESKY_DECOMP(I_mat);
  gsl_linalg_cholesky_invert(I_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);
}

/* =========================================================
   LEGACY FUNCTION — ridgeReg (.C interface, 32-bit)

   Optimal for single-threaded execution on any platform.
   Best default on macOS + vecLib for typical dataset sizes.
   Uses Y-row permutation with contiguous row copies into
   Y_rand, then one large dgemm(T, Y_rand) per permutation.
   This maximises BLAS efficiency in the single-threaded case.
   ========================================================= */
void ridgeReg(
  double *X_vec, double *Y_vec,
  int *n_pt, int *p_pt, int *m_pt,
  double *lambda_pt, double *nrand_pt,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int *rng_method_pt)
{
  gsl_matrix *X, *Y, *I, *T, *beta, *Y_rand, *beta_rand, *aver, *aver_sq, *zscore, *pvalue;
  int n = *n_pt, p = *p_pt, m = *m_pt, nrand = (int)*nrand_pt;
  double lambda = *lambda_pt;
  int rng_method = *rng_method_pt;
  int *array_index, i, j;

  X       = RVectorObject_to_gsl_matrix(X_vec,     n, p);
  Y       = RVectorObject_to_gsl_matrix(Y_vec,     n, m);
  beta    = RVectorObject_to_gsl_matrix(beta_vec,  p, m);
  aver_sq = RVectorObject_to_gsl_matrix(se_vec,    p, m);
  zscore  = RVectorObject_to_gsl_matrix(zscore_vec,p, m);
  pvalue  = RVectorObject_to_gsl_matrix(pvalue_vec,p, m);

  I         = gsl_matrix_alloc(p, p);
  T         = gsl_matrix_alloc(p, n);
  Y_rand    = gsl_matrix_alloc(n, m);
  beta_rand = gsl_matrix_alloc(p, m);
  aver      = gsl_matrix_alloc(p, m);
  array_index = (int*)malloc(n * sizeof(int));

  /* Projection: T = (X'X + lambda*I)^-1 * X' */
  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);
  CHOLESKY_DECOMP(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);

  /* Observed beta */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

  for (i = 0; i < n; i++) array_index[i] = i;
  gsl_matrix_set_zero(aver);
  gsl_matrix_set_zero(aver_sq);
  gsl_matrix_set_zero(pvalue);

  gsl_rng *gsl_rng_obj = NULL;
  if (rng_method == RNG_GSL) {
    gsl_rng_obj = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gsl_rng_obj, 0);
  } else {
    srand(0);
  }

  for (i = 0; i < nrand; i++) {
    /* Shuffle row indices */
    if (rng_method == RNG_GSL) shuffle_gsl(gsl_rng_obj, array_index, n);
    else                        shuffle_srand(array_index, n);

    /* Contiguous row copy: Y_rand[j,:] = Y[array_index[j],:]
       Each row is m=1000 contiguous doubles — very cache-friendly. */
    for (j = 0; j < n; j++) {
      gsl_vector_const_view t = gsl_matrix_const_row(Y, array_index[j]);
      gsl_matrix_set_row(Y_rand, j, &t.vector);
    }

    /* One large dgemm per permutation — maximally BLAS-efficient. */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);

    for (j = 0; j < (int)(pvalue->size1 * pvalue->size2); j++) {
      if (fabs(beta_rand->data[j]) >= fabs(beta->data[j])) pvalue->data[j]++;
    }
    gsl_matrix_add(aver, beta_rand);
    gsl_matrix_mul_elements(beta_rand, beta_rand);
    gsl_matrix_add(aver_sq, beta_rand);
  }

  gsl_matrix_scale(aver,    1.0 / nrand);
  gsl_matrix_scale(aver_sq, 1.0 / nrand);
  gsl_matrix_add_constant(pvalue, 1.0);
  gsl_matrix_scale(pvalue, 1.0 / (nrand + 1));

  gsl_matrix_memcpy(zscore, beta);
  gsl_matrix_sub(zscore, aver);
  gsl_matrix_mul_elements(aver, aver);
  gsl_matrix_sub(aver_sq, aver);

  for (i = 0; i < (int)(aver_sq->size1 * aver_sq->size2); i++)
    aver_sq->data[i] = sqrt(aver_sq->data[i] > 0 ? aver_sq->data[i] : 0);

  gsl_matrix_div_elements(zscore, aver_sq);

  if (gsl_rng_obj) gsl_rng_free(gsl_rng_obj);
  gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_rand);
  gsl_matrix_free(beta_rand); gsl_matrix_free(aver);
  free(array_index);
  gsl_matrix_partial_free(beta); gsl_matrix_partial_free(aver_sq);
  gsl_matrix_partial_free(zscore); gsl_matrix_partial_free(pvalue);
}


/* =========================================================
   YROW CORE (.Call, 64-bit)

   Single-threaded path (num_threads == 1):
     Mirrors ridgeReg's optimal design.
     One-time transpose of Yt(m×n) → Y_real(n×m) so that
     row copies into Y_rand are contiguous (cache-friendly).
     One large dgemm(T, Y_rand) per permutation.
     No batching overhead.

   Multi-threaded path (num_threads > 1):
     Sample-strip + perm-batch approach. Each thread owns a
     fixed-size working set. Parallelises over sample strips.
     OpenMP disabled on macOS+vecLib via configure.
   ========================================================= */
void ridgeRegFast_core(
  double *X_ptr, double *Y_ptr,
  size_t n, size_t p, size_t m,
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int num_threads, int rng_method)
{
  /* Defer omp_set_num_threads() to just before the parallel region.
     Calling it for the single-threaded path (num_threads == 1)
     triggers libgomp CPU-affinity pinning, which forces OpenBLAS's
     internal pthreads onto one core → catastrophic slowdown. */

  gsl_matrix *Xt   = RVectorObject_to_gsl_matrix(X_ptr,    p, n);
  gsl_matrix *Yt   = RVectorObject_to_gsl_matrix(Y_ptr,    m, n);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T     = gsl_matrix_alloc(p, n);

  compute_projection(Xt, I_mat, T, p, lambda);

  /* Observed beta = T * Y_real = T * Yt' */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  /* Build forward permutation table (shared between st and mt paths) */
  int *perm_table = build_perm_table(nrand, n, rng_method);

  /* Initialise output accumulators */
  R_xlen_t total = (R_xlen_t)p * (R_xlen_t)m;
  for (R_xlen_t i = 0; i < total; i++) {
    pvalue_vec[i] = 0.0; se_vec[i] = 0.0; zscore_vec[i] = 0.0;
  }

  /* -------------------------------------------------------
     SINGLE-THREADED PATH (num_threads == 1)
     -------------------------------------------------------
     Transpose Yt(m×n) → Y_real(n×m) once.
     Row j of Y_real (m contiguous doubles) = sample values
     for gene j — identical layout to ridgeReg's Y matrix.
     Per permutation: cache-friendly row copies into Y_rand,
     then one maximally-sized dgemm(T, Y_rand).
     No batching, no strip overhead.
     ------------------------------------------------------- */
  if (num_threads == 1) {
    /* One-time transpose: Yt(m×n) → Y_real(n×m).
       Cost: O(n×m) — amortised over nrand permutations. */
    gsl_matrix *Y_real = gsl_matrix_alloc(n, m);
    for (size_t gene = 0; gene < n; gene++)
      for (size_t samp = 0; samp < m; samp++)
        gsl_matrix_set(Y_real, gene, samp, gsl_matrix_get(Yt, samp, gene));

    gsl_matrix *Y_rand    = gsl_matrix_alloc(n, m);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

    for (int i_rand = 0; i_rand < nrand; i_rand++) {
      int *perm = &perm_table[(size_t)i_rand * n];

      /* Contiguous row copy: Y_rand[j,:] = Y_real[perm[j],:]
         Each row is m doubles — cache-line friendly. */
      for (size_t j = 0; j < n; j++) {
        gsl_vector_const_view src = gsl_matrix_const_row(Y_real, (size_t)perm[j]);
        gsl_matrix_set_row(Y_rand, j, &src.vector);
      }

      /* One large dgemm per permutation — optimal BLAS efficiency */
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);

      accumulate_stats_serial(beta_vec, se_vec, zscore_vec, pvalue_vec,
                              beta_rand, p, m);
    }

    gsl_matrix_free(Y_real);
    gsl_matrix_free(Y_rand);
    gsl_matrix_free(beta_rand);

  } else {
    /* -------------------------------------------------------
       MULTI-THREADED PATH (num_threads > 1)
       -------------------------------------------------------
       Parallelises over sample strips. Each thread maintains
       a fixed-size working set (SAMPLE_STRIP_SIZE × PERM_BATCH_SIZE)
       regardless of total m. Permutation batching amortises
       the Y_block construction overhead across threads.
       ------------------------------------------------------- */
#ifdef _OPENMP
    if (num_threads > 0) omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
      gsl_matrix *Y_block    = gsl_matrix_alloc(n, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
      gsl_matrix *Beta_block = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);

      #pragma omp for schedule(dynamic)
      for (size_t samp_start = 0; samp_start < m; samp_start += SAMPLE_STRIP_SIZE) {
        size_t cur_samp = (samp_start + SAMPLE_STRIP_SIZE > m) ? (m - samp_start) : SAMPLE_STRIP_SIZE;

        for (int b_start = 0; b_start < nrand; b_start += PERM_BATCH_SIZE) {
          int cur_perm = (b_start + PERM_BATCH_SIZE > nrand) ? (nrand - b_start) : PERM_BATCH_SIZE;

          for (int i_perm = 0; i_perm < cur_perm; i_perm++) {
            int *p_idx = &perm_table[(size_t)(b_start + i_perm) * n];
            for (size_t s_local = 0; s_local < cur_samp; s_local++) {
              double *src_row = gsl_matrix_ptr(Yt, samp_start + s_local, 0);
              size_t  dest_col = (size_t)i_perm * cur_samp + s_local;
              for (size_t gene = 0; gene < n; gene++)
                gsl_matrix_set(Y_block, gene, dest_col, src_row[p_idx[gene]]);
            }
          }

          gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block,    0, 0, n, cur_perm * cur_samp);
          gsl_matrix_view B_sub = gsl_matrix_submatrix(Beta_block, 0, 0, p, cur_perm * cur_samp);
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);

          for (int i_perm = 0; i_perm < cur_perm; i_perm++) {
            for (size_t r = 0; r < p; r++) {
              for (size_t s_local = 0; s_local < cur_samp; s_local++) {
                R_xlen_t idx   = (R_xlen_t)r * m + (samp_start + s_local);
                double   b_rnd = gsl_matrix_get(&B_sub.matrix, r, (size_t)i_perm * cur_samp + s_local);
                double   b_obs = beta_vec[idx];
                #pragma omp atomic
                pvalue_vec[idx] += (fabs(b_rnd) >= fabs(b_obs)) ? 1.0 : 0.0;
                #pragma omp atomic
                zscore_vec[idx] += b_rnd;
                #pragma omp atomic
                se_vec[idx]     += b_rnd * b_rnd;
              }
            }
          }
        }
      }
      gsl_matrix_free(Y_block);
      gsl_matrix_free(Beta_block);
    }
#endif
  }

  free(perm_table);
  finalize_stats(beta_vec, se_vec, zscore_vec, pvalue_vec, total, nrand);

  gsl_matrix_free(I_mat); gsl_matrix_free(T);
  gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt);
  gsl_matrix_partial_free(beta);
}


/* =========================================================
   TCOL CORE (.Call, 64-bit)

   Single-threaded path (num_threads == 1):
     Per permutation: construct T_perm(p×n) by permuting
     columns of T — small p=500 doubles per column, fully
     cache-resident. One large dgemm(T_perm, Yt^T) against
     the full Y. No sample strips, no batching overhead.
     Substantially faster than the mt strip approach for
     single-threaded execution.

   Multi-threaded path (num_threads > 1):
     Parallelises over permutations. Each thread constructs
     its own T_perm and processes all sample strips for that
     permutation. Atomic accumulation into shared output.
   ========================================================= */
void ridgeRegFastTcol_core(
  double *X_ptr, double *Y_ptr,
  size_t n, size_t p, size_t m,
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int num_threads, int rng_method)
{
  /* Defer omp_set_num_threads() to just before the parallel region.
     See Yrow core comment for rationale. */

  gsl_matrix *Xt   = RVectorObject_to_gsl_matrix(X_ptr,    p, n);
  gsl_matrix *Yt   = RVectorObject_to_gsl_matrix(Y_ptr,    m, n);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T     = gsl_matrix_alloc(p, n);

  compute_projection(Xt, I_mat, T, p, lambda);

  /* Observed beta = T * Yt' */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  /* Build INVERSE permutation table for T-column approach.
     Forward shuffle → inv_perm[fwd[j]] = j.
     T_perm[:,col] = T[:,inv_perm[col]] mathematically equivalent
     to permuting rows of Y. */
  int *fwd_table = build_perm_table(nrand, n, rng_method);
  int *inv_perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  for (int i_rand = 0; i_rand < nrand; i_rand++) {
    int *fwd = &fwd_table[(size_t)i_rand * n];
    int *inv = &inv_perm_table[(size_t)i_rand * n];
    for (size_t j = 0; j < n; j++) inv[fwd[j]] = (int)j;
  }
  free(fwd_table);

  /* Initialise output accumulators */
  R_xlen_t total = (R_xlen_t)p * (R_xlen_t)m;
  for (R_xlen_t i = 0; i < total; i++) {
    pvalue_vec[i] = 0.0; se_vec[i] = 0.0; zscore_vec[i] = 0.0;
  }

  /* -------------------------------------------------------
     SINGLE-THREADED PATH (num_threads == 1)
     -------------------------------------------------------
     For each permutation:
       1. Build T_perm(p×n) by permuting columns of T.
          Column width = p = 500 doubles = 4 KB — fits in L1.
          n=16325 columns → T_perm ≈ 65 MB, accessed once
          sequentially → cache-friendly.
       2. One large dgemm(T_perm, Yt^T) → beta_rand(p×m).
          Shape: (500×16325) × (16325×1000) — maximally
          BLAS-efficient, same as ridgeReg's dgemm shape.
     No sample strips, no atomics, no batch overhead.
     ------------------------------------------------------- */
  if (num_threads == 1) {
    gsl_matrix *T_perm    = gsl_matrix_alloc(p, n);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

    for (int i_rand = 0; i_rand < nrand; i_rand++) {
      int *inv_perm = &inv_perm_table[(size_t)i_rand * n];

      /* Permute columns of T → T_perm.
         Each column: p=500 doubles, stride-1 within column.
         Sequential column iteration → predictable prefetch. */
      for (size_t col = 0; col < n; col++) {
        size_t src_col = (size_t)inv_perm[col];
        for (size_t row = 0; row < p; row++)
          gsl_matrix_set(T_perm, row, col,
                         gsl_matrix_get(T, row, src_col));
      }

      /* One large dgemm against the full Y — optimal shape */
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T_perm, Yt, 0.0, beta_rand);

      accumulate_stats_serial(beta_vec, se_vec, zscore_vec, pvalue_vec,
                              beta_rand, p, m);
    }

    gsl_matrix_free(T_perm);
    gsl_matrix_free(beta_rand);

  } else {
    /* -------------------------------------------------------
       MULTI-THREADED PATH (num_threads > 1)
       -------------------------------------------------------
       Parallelises over permutations. Each thread owns a
       private T_perm and beta_strip. Atomic accumulation
       into shared pvalue/zscore/se output arrays.
       ------------------------------------------------------- */
#ifdef _OPENMP
    if (num_threads > 0) omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
      gsl_matrix *T_perm     = gsl_matrix_alloc(p, n);
      gsl_matrix *beta_strip = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE);

      #pragma omp for schedule(dynamic)
      for (int i_rand = 0; i_rand < nrand; i_rand++) {
        int *inv_perm = &inv_perm_table[(size_t)i_rand * n];

        /* Permute columns of T */
        for (size_t col = 0; col < n; col++) {
          size_t src_col = (size_t)inv_perm[col];
          for (size_t row = 0; row < p; row++)
            gsl_matrix_set(T_perm, row, col,
                           gsl_matrix_get(T, row, src_col));
        }

        /* Process sample strips */
        for (size_t samp_start = 0; samp_start < m; samp_start += SAMPLE_STRIP_SIZE) {
          size_t cur_samp = (samp_start + SAMPLE_STRIP_SIZE > m)
                             ? (m - samp_start) : SAMPLE_STRIP_SIZE;

          gsl_matrix_view Y_sub = gsl_matrix_submatrix(Yt, samp_start, 0, cur_samp, n);
          gsl_matrix_view B_sub = gsl_matrix_submatrix(beta_strip, 0, 0, p, cur_samp);
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T_perm, &Y_sub.matrix, 0.0, &B_sub.matrix);

          for (size_t r = 0; r < p; r++) {
            for (size_t s_local = 0; s_local < cur_samp; s_local++) {
              R_xlen_t idx   = (R_xlen_t)r * m + (samp_start + s_local);
              double   b_rnd = gsl_matrix_get(&B_sub.matrix, r, s_local);
              double   b_obs = beta_vec[idx];
              #pragma omp atomic
              pvalue_vec[idx] += (fabs(b_rnd) >= fabs(b_obs)) ? 1.0 : 0.0;
              #pragma omp atomic
              zscore_vec[idx] += b_rnd;
              #pragma omp atomic
              se_vec[idx]     += b_rnd * b_rnd;
            }
          }
        }
      }
      gsl_matrix_free(T_perm);
      gsl_matrix_free(beta_strip);
    }
#endif
  }

  free(inv_perm_table);
  finalize_stats(beta_vec, se_vec, zscore_vec, pvalue_vec, total, nrand);

  gsl_matrix_free(I_mat); gsl_matrix_free(T);
  gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt);
  gsl_matrix_partial_free(beta);
}


/* =========================================================
   .Call INTERFACES
   ========================================================= */

typedef void (*ridge_core_fn)(double*, double*, size_t, size_t, size_t,
                               double, int, double*, double*, double*, double*,
                               int, int);

static SEXP ridge_interface_common(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp,
                                    SEXP nrand_sexp, SEXP ncores_sexp,
                                    SEXP rng_method_sexp, ridge_core_fn core_fn)
{
  SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
  size_t n = (size_t)INTEGER(x_dim)[0];
  size_t p = (size_t)INTEGER(x_dim)[1];
  SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
  size_t m = (size_t)INTEGER(y_dim)[1];
  R_xlen_t total = (R_xlen_t)p * (R_xlen_t)m;

  SEXP beta_s   = PROTECT(allocVector(REALSXP, total));
  SEXP se_s     = PROTECT(allocVector(REALSXP, total));
  SEXP zscore_s = PROTECT(allocVector(REALSXP, total));
  SEXP pvalue_s = PROTECT(allocVector(REALSXP, total));

  core_fn(REAL(X_sexp), REAL(Y_sexp), n, p, m,
          asReal(lambda_sexp), asInteger(nrand_sexp),
          REAL(beta_s), REAL(se_s), REAL(zscore_s), REAL(pvalue_s),
          asInteger(ncores_sexp), asInteger(rng_method_sexp));

  SEXP res   = PROTECT(allocVector(VECSXP, 4));
  SEXP names = PROTECT(allocVector(STRSXP, 4));
  SET_VECTOR_ELT(res, 0, beta_s);   SET_STRING_ELT(names, 0, mkChar("beta"));
  SET_VECTOR_ELT(res, 1, se_s);     SET_STRING_ELT(names, 1, mkChar("se"));
  SET_VECTOR_ELT(res, 2, zscore_s); SET_STRING_ELT(names, 2, mkChar("zscore"));
  SET_VECTOR_ELT(res, 3, pvalue_s); SET_STRING_ELT(names, 3, mkChar("pvalue"));
  setAttrib(res, R_NamesSymbol, names);
  UNPROTECT(6);
  return res;
}

SEXP ridgeRegFast_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp,
                             SEXP nrand_sexp, SEXP ncores_sexp, SEXP rng_method_sexp)
{
  return ridge_interface_common(X_sexp, Y_sexp, lambda_sexp,
                                 nrand_sexp, ncores_sexp, rng_method_sexp,
                                 ridgeRegFast_core);
}

SEXP ridgeRegFastTcol_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp,
                                 SEXP nrand_sexp, SEXP ncores_sexp, SEXP rng_method_sexp)
{
  return ridge_interface_common(X_sexp, Y_sexp, lambda_sexp,
                                 nrand_sexp, ncores_sexp, rng_method_sexp,
                                 ridgeRegFastTcol_core);
}


/* =========================================================
   GENERATE PERMUTATION TABLE (.Call)
   Returns integer matrix (nrand x n) of forward permutation
   indices. Uses same RNG sequence as ridgeReg/ridgeRegFast.
   Enables the naive R variant to match C RNG exactly.
   ========================================================= */
SEXP generate_perm_table(SEXP n_sexp, SEXP nrand_sexp, SEXP rng_method_sexp)
{
  int n      = asInteger(n_sexp);
  int nrand  = asInteger(nrand_sexp);
  int rng_method = asInteger(rng_method_sexp);

  SEXP result = PROTECT(allocMatrix(INTSXP, nrand, n));
  int *result_ptr = INTEGER(result);

  int *perm_table = build_perm_table(nrand, (size_t)n, rng_method);

  /* Store as column-major R matrix (nrand×n): element [i,j] at j*nrand+i */
  for (int i = 0; i < nrand; i++)
    for (int j = 0; j < n; j++)
      result_ptr[(size_t)j * nrand + i] = perm_table[(size_t)i * n + j];

  free(perm_table);
  UNPROTECT(1);
  return result;
}


/* =========================================================
   REGISTRATION
   ========================================================= */
static const R_CMethodDef cMethods[] = {
  {"ridgeReg", (DL_FUNC) &ridgeReg, 12},
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
