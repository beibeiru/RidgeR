#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <string.h> 
#include <time.h>   
#include <math.h>   
#include <stdlib.h> 

#ifdef _OPENMP
  #include <omp.h>
#endif

// --- DYNAMIC TUNING CONSTANTS ---
// We target a buffer size that fits comfortably in L2 Cache per thread.
// 2 MB = ~262,144 doubles.
// This prevents memory bandwidth saturation on Apple Silicon.
#define TARGET_CACHE_ELEMENTS 262144 

/* =========================================================
   HELPER FUNCTIONS
   ========================================================= */

void disable_blas_threading() {
    #ifdef _WIN32
        _putenv("OPENBLAS_NUM_THREADS=1");
        _putenv("MKL_NUM_THREADS=1");
        _putenv("VECLIB_MAXIMUM_THREADS=1");
    #else
        setenv("OPENBLAS_NUM_THREADS", "1", 1);
        setenv("MKL_NUM_THREADS", "1", 1);
        setenv("VECLIB_MAXIMUM_THREADS", "1", 1);
    #endif
}

gsl_matrix *RVectorObject_to_gsl_matrix(double *vec, size_t nr, size_t nc) {
  gsl_block *b = (gsl_block*)malloc(sizeof(gsl_block));
  gsl_matrix *r = (gsl_matrix*)malloc(sizeof(gsl_matrix));
  r->size1 = nr; r->tda = r->size2 = nc; r->owner = 1; 
  b->data = r->data = vec; r->block = b; b->size = r->size1 * r->size2;
  return r;
}

void gsl_matrix_partial_free(gsl_matrix *x) {
  if(x) { if(x->block) free(x->block); free(x); }
}

void shuffle(int array[], const int n) {
  int i, j; double t;
  for (i = 0; i < n-1; i++) {
    j = i + rand() / (RAND_MAX / (n - i) + 1);
    t = array[j]; array[j] = array[i]; array[i] = t;
  }
}

/* =========================================================
   OLD VERSION (Legacy)
   ========================================================= */
void ridgeReg(
  double *X_vec, double *Y_vec,
  int *n_pt, int *p_pt, int *m_pt,
  double *lambda_pt, double *nrand_pt,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec
) {
  gsl_matrix *X = RVectorObject_to_gsl_matrix(X_vec, *n_pt, *p_pt);
  gsl_matrix *Y = RVectorObject_to_gsl_matrix(Y_vec, *n_pt, *m_pt);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, *p_pt, *m_pt);
  gsl_matrix *aver_sq = RVectorObject_to_gsl_matrix(se_vec, *p_pt, *m_pt);
  gsl_matrix *zscore = RVectorObject_to_gsl_matrix(zscore_vec, *p_pt, *m_pt);
  gsl_matrix *pvalue = RVectorObject_to_gsl_matrix(pvalue_vec, *p_pt, *m_pt);

  int n = *n_pt, p = *p_pt, m = *m_pt, nrand = (int)*nrand_pt;
  double lambda = *lambda_pt;

  gsl_matrix *I = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);
  gsl_matrix *Y_rand = gsl_matrix_alloc(n, m);
  gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
  gsl_matrix *aver = gsl_matrix_alloc(p, m);
  int *array_index = (int*)malloc(n*sizeof(int));

  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

  srand(0);
  for(int i=0;i<n;i++) array_index[i] = i;
  gsl_matrix_set_zero(aver); gsl_matrix_set_zero(aver_sq); gsl_matrix_set_zero(pvalue);

  for(int i=0;i<nrand;i++) {
    shuffle(array_index, n);
    for(int j=0;j<n;j++){
      gsl_vector_const_view t = gsl_matrix_const_row(Y, array_index[j]);
      gsl_matrix_set_row(Y_rand, j, &t.vector);
    }
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);
    for(int j=0; j< p*m; j++) {
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
  for(int i=0;i< p*m; i++) aver_sq->data[i] = sqrt(aver_sq->data[i]);
  gsl_matrix_div_elements(zscore, aver_sq);

  gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_rand);
  gsl_matrix_free(beta_rand); gsl_matrix_free(aver); free(array_index);
  gsl_matrix_partial_free(beta); gsl_matrix_partial_free(aver_sq);
  gsl_matrix_partial_free(zscore); gsl_matrix_partial_free(pvalue);
}

/* =========================================================
   NEW OPTIMIZED LOGIC (Dynamic Cache Tuning)
   ========================================================= */

void ridgeRegFast_core(
  double *X_ptr, double *Y_ptr,
  size_t n, size_t p, size_t m,
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int num_threads
)
{
  if (num_threads > 0) {
       #ifdef _OPENMP
           omp_set_num_threads(num_threads);
           disable_blas_threading();
       #endif
  }

  // --- DYNAMIC TUNING ---
  // We calculate block sizes to stay within L2 cache limits.
  // Formula: n * strip_size * batch_size <= TARGET_CACHE_ELEMENTS
  
  // 1. Determine optimal Strip Size (Samples per block)
  // We want at least 16 samples for vectorization, max 128 for responsiveness.
  size_t sample_strip_size = 64; 
  
  // 2. Determine optimal Batch Size (Permutations per block) based on remaining memory
  size_t elements_per_perm_strip = n * sample_strip_size;
  size_t perm_batch_size = 1;
  
  if (elements_per_perm_strip < TARGET_CACHE_ELEMENTS) {
      perm_batch_size = TARGET_CACHE_ELEMENTS / elements_per_perm_strip;
  }
  
  // Safety Clamps
  if (perm_batch_size < 1) perm_batch_size = 1;
  if (perm_batch_size > 32) perm_batch_size = 32; // Diminishing returns after 32

  // --- SETUP ---
  gsl_matrix *X_gsl = RVectorObject_to_gsl_matrix(X_ptr, n, p);
  gsl_matrix *Y_gsl = RVectorObject_to_gsl_matrix(Y_ptr, n, m);
  gsl_matrix *beta  = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);

  gsl_matrix_set_identity(I_mat);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X_gsl, lambda, I_mat);
  gsl_linalg_cholesky_decomp(I_mat);
  gsl_linalg_cholesky_invert(I_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I_mat, X_gsl, 0.0, T);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_gsl, 0.0, beta);

  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  for(size_t k=0; k<n; k++) temp_idx[k] = (int)k;

  srand(0);
  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, (int)n);
    memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
  }
  free(temp_idx);

  R_xlen_t total_elements = (R_xlen_t)p * (R_xlen_t)m;
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = 0.0; se_vec[i] = 0.0; zscore_vec[i] = 0.0; 
  }

  // --- PARALLEL LOOP ---
  #pragma omp parallel
  {
    // Buffers allocated per thread
    gsl_matrix *Y_block = gsl_matrix_alloc(n, sample_strip_size * perm_batch_size);
    gsl_matrix *Beta_block = gsl_matrix_alloc(p, sample_strip_size * perm_batch_size);
    
    #pragma omp for schedule(dynamic)
    for(size_t samp_start = 0; samp_start < m; samp_start += sample_strip_size) 
    {
      size_t current_samples = (samp_start + sample_strip_size > m) ? (m - samp_start) : sample_strip_size;
      
      for(int b_start = 0; b_start < nrand; b_start += (int)perm_batch_size)
      {
        int current_perms = (b_start + (int)perm_batch_size > nrand) ? (nrand - b_start) : (int)perm_batch_size;
        size_t total_cols = current_perms * current_samples;

        // A. Construct Block (Fast Memcpy)
        for(int i_perm = 0; i_perm < current_perms; i_perm++) {
           int *p_idx = &perm_table[(size_t)(b_start + i_perm) * n];
           size_t dest_col_offset = i_perm * current_samples;
           
           for(size_t row = 0; row < n; row++) {
             size_t src_row_idx = p_idx[row];
             double *src_ptr = gsl_matrix_ptr(Y_gsl, src_row_idx, samp_start);
             double *dest_ptr = gsl_matrix_ptr(Y_block, row, dest_col_offset);
             memcpy(dest_ptr, src_ptr, current_samples * sizeof(double));
           }
        }

        // B. Matrix Multiply
        gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block, 0, 0, n, total_cols);
        gsl_matrix_view B_sub = gsl_matrix_submatrix(Beta_block, 0, 0, p, total_cols);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);

        // C. Update Stats
        for(int i_perm = 0; i_perm < current_perms; i_perm++) {
          size_t col_offset = i_perm * current_samples;
          
          for(size_t r = 0; r < p; r++) {
            double *beta_computed = gsl_matrix_ptr(&B_sub.matrix, r, col_offset);
            
            for(size_t s_local = 0; s_local < current_samples; s_local++) {
              R_xlen_t global_idx = (R_xlen_t)r * m + (samp_start + s_local);
              double b_rnd = beta_computed[s_local];
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
  gsl_matrix_partial_free(X_gsl); gsl_matrix_partial_free(Y_gsl); gsl_matrix_partial_free(beta);
}

// .Call Interface
SEXP ridgeRegFast_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp, SEXP ncores_sexp) 
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

    ridgeRegFast_core(
        REAL(X_sexp), REAL(Y_sexp), n, p, m,
        asReal(lambda_sexp), asInteger(nrand_sexp),
        REAL(beta_sexp), REAL(se_sexp), REAL(zscore_sexp), REAL(pvalue_sexp),
        asInteger(ncores_sexp)
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

// Function Registration
static const R_CMethodDef cMethods[] = {
    {"ridgeReg", (DL_FUNC) &ridgeReg, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef callMethods[] = {
    {"ridgeRegFast_interface", (DL_FUNC) &ridgeRegFast_interface, 5},
    {NULL, NULL, 0}
};

void R_init_SecAct(DllInfo *dll) { 
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
