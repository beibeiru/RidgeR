#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

// --- TUNING PARAMETERS ---
#define SAMPLE_STRIP_SIZE 64 
#define PERM_BATCH_SIZE 64

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

void shuffle(int array[], const int n)
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

/* =========================================================
   VERSION 1: ORIGINAL (Single-threaded, Y row permutation)
   Used by: SecAct.inference.gsl.old
   ========================================================= */
void ridgeReg(
  double *X_vec, double *Y_vec,
  int *n_pt, int *p_pt, int *m_pt,
  double *lambda_pt, double *nrand_pt,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec
)
{
  gsl_matrix *X, *Y, *I, *T, *beta, *Y_rand, *beta_rand, *aver, *aver_sq, *zscore, *pvalue;
  int n = *n_pt, p = *p_pt, m = *m_pt, nrand = (int)*nrand_pt;
  double lambda = *lambda_pt;
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

  // Compute T = (X'X + lambda*I)^-1 * X'
  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);
  
  // Compute observed beta = T * Y
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

  srand(0);
  for(i=0;i<n;i++) array_index[i] = i;

  gsl_matrix_set_zero(aver);
  gsl_matrix_set_zero(aver_sq);
  gsl_matrix_set_zero(pvalue);

  // Permutation test: permute rows of Y
  for(i=0;i<nrand;i++)
  {
    shuffle(array_index, n);

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
   VERSION 1.5: SINGLE THREADED T-PERM (Cache Optimized)
   Used by: SecAct.inference.gsl.old2
   ========================================================= */
void ridgeRegTperm_single(
  double *X_vec, double *Y_vec,
  int *n_pt, int *p_pt, int *m_pt,
  double *lambda_pt, double *nrand_pt,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec
)
{
  int n = *n_pt;
  int p = *p_pt;
  int m = *m_pt;
  int nrand = (int)*nrand_pt;
  double lambda = *lambda_pt;

  // Matrix views (using the "Transposed view" trick for correct math)
  // Xt: (p x n) view of column-major X
  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_vec, p, n);
  // Yt: (m x n) view of column-major Y
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_vec, m, n);
  // beta: (p x m)
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);

  // T = (X'X + lambda*I)^-1 * X'
  gsl_matrix_set_identity(I_mat);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
  gsl_linalg_cholesky_decomp(I_mat);
  gsl_linalg_cholesky_invert(I_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);

  // Observed beta = T * Y
  // Math: beta = T @ Y
  // GSL:  beta = T @ Yt^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  // --- PREPARE FOR CACHE-EFFICIENT PERMUTATION ---
  // Create T^T (n x p). GSL stores Row-Major.
  gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
  gsl_matrix_transpose_memcpy(Tt_orig, T);

  gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);
  gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

  int *p_idx = (int*)malloc(n * sizeof(int));
  for(int k=0; k<n; k++) p_idx[k] = k;

  // Helpers for memcpy
  double *src_base = Tt_orig->data;
  double *dst_base = Tt_perm->data;
  size_t row_size_bytes = p * sizeof(double);

  // Initialize Statistics
  size_t total_elements = (size_t)p * m;
  for(size_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = 0.0;
    se_vec[i] = 0.0;
    zscore_vec[i] = 0.0;
  }

  srand(0);

  for(int i_rand=0; i_rand<nrand; i_rand++)
  {
    shuffle(p_idx, n);

    // Optimized Permutation: Row-wise memcpy on T^T
    for(int r=0; r<n; r++) {
       memcpy(dst_base + (r * p), 
              src_base + (p_idx[r] * p), 
              row_size_bytes);
    }

    // beta_rand = T_perm @ Y
    //           = (Tt_perm)^T @ Yt^T
    gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Tt_perm, Yt, 0.0, beta_rand);

    // Accumulate
    for(size_t k=0; k<total_elements; k++) {
      double b_rnd = beta_rand->data[k];
      double b_obs = beta_vec[k];

      if(fabs(b_rnd) >= fabs(b_obs)) pvalue_vec[k] += 1.0;
      zscore_vec[k] += b_rnd;
      se_vec[k] += b_rnd * b_rnd;
    }
  }

  // Finalize Statistics
  double inv_nrand = 1.0 / nrand;
  for(size_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean_rand = zscore_vec[i] * inv_nrand;
    double var_rand = (se_vec[i] * inv_nrand) - (mean_rand * mean_rand);
    double se_rand = sqrt(var_rand > 0 ? var_rand : 0);
    se_vec[i] = se_rand;
    zscore_vec[i] = (se_rand > 1e-12) ? (beta_vec[i] - mean_rand) / se_rand : 0.0;
  }

  free(p_idx);
  gsl_matrix_free(Tt_orig);
  gsl_matrix_free(Tt_perm);
  gsl_matrix_free(beta_rand);
  gsl_matrix_free(I_mat);
  gsl_matrix_free(T);
  gsl_matrix_partial_free(Xt);
  gsl_matrix_partial_free(Yt);
  gsl_matrix_partial_free(beta);
}

/* =========================================================
   VERSION 2: OPTIMIZED (Multi-threaded, Y row permutation with blocking)
   Used by: SecAct.inference.gsl.new
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
    #endif
  }

  // R stores X as (n x p) Column-Major.
  // GSL reads Row-Major. 
  // Mapping X_ptr as (p x n) gives us t(X) directly.
  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n);
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);

  // T = (X'X + lambda*I)^-1 * X'
  gsl_matrix_set_identity(I_mat);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
  gsl_linalg_cholesky_decomp(I_mat);
  gsl_linalg_cholesky_invert(I_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);

  // beta = T * Y (using transposed views: T @ Yt^T)
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  // Pre-generate permutation table
  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  for(size_t k=0; k<n; k++) temp_idx[k] = (int)k;

  srand(0);
  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, (int)n);
    memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
  }
  free(temp_idx);

  // Initialize outputs
  R_xlen_t total_elements = (R_xlen_t)p * (R_xlen_t)m;
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = 0.0; se_vec[i] = 0.0; zscore_vec[i] = 0.0; 
  }

  // Parallel loop with blocking over samples
  #pragma omp parallel
  {
    gsl_matrix *Y_block = gsl_matrix_alloc(n, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    gsl_matrix *Beta_block = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    
    #pragma omp for schedule(dynamic)
    for(size_t samp_start = 0; samp_start < m; samp_start += SAMPLE_STRIP_SIZE) 
    {
      size_t current_samples = (samp_start + SAMPLE_STRIP_SIZE > m) ? (m - samp_start) : SAMPLE_STRIP_SIZE;
      
      for(int b_start = 0; b_start < nrand; b_start += PERM_BATCH_SIZE)
      {
        int current_perms = (b_start + PERM_BATCH_SIZE > nrand) ? (nrand - b_start) : PERM_BATCH_SIZE;

        // Construct Y_block with permuted rows
        for(int i_perm = 0; i_perm < current_perms; i_perm++) {
          int *p_idx = &perm_table[(size_t)(b_start + i_perm) * n];
          
          for(size_t s_local = 0; s_local < current_samples; s_local++) {
            size_t global_samp_idx = samp_start + s_local;
            double *src_row = gsl_matrix_ptr(Yt, global_samp_idx, 0);
            size_t dest_col = (i_perm * current_samples) + s_local;
            
            for(size_t gene = 0; gene < n; gene++) {
              double val = src_row[ p_idx[gene] ];
              gsl_matrix_set(Y_block, gene, dest_col, val);
            }
          }
        }

        gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block, 0, 0, n, current_perms * current_samples);
        gsl_matrix_view B_sub = gsl_matrix_submatrix(Beta_block, 0, 0, p, current_perms * current_samples);
        
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);

        // Accumulate statistics
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
    gsl_matrix_free(Y_block); 
    gsl_matrix_free(Beta_block);
  }

  free(perm_table);

  // Finalize statistics
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

  gsl_matrix_free(I_mat); 
  gsl_matrix_free(T);
  gsl_matrix_partial_free(Xt); 
  gsl_matrix_partial_free(Yt); 
  gsl_matrix_partial_free(beta);
}

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

/* =========================================================
   VERSION 3: T COLUMN PERMUTATION (Multi-threaded)
   Used by: SecAct.inference.gsl.new2
   ========================================================= */

void ridgeRegTperm_core(
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
    #endif
  }

  // Matrix views
  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n);  // t(X): p x n
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n);  // t(Y): m x n (GSL view of Col-Major R matrix)
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);

  // T = (X'X + lambda*I)^-1 * X'
  gsl_matrix_set_identity(I_mat);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
  gsl_linalg_cholesky_decomp(I_mat);
  gsl_linalg_cholesky_invert(I_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);

  // Observed beta = T @ Y (= T @ Yt^T in our view)
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  // --- PREPARE FOR CACHE-EFFICIENT PERMUTATION ---
  // Create T^T (n x p). GSL stores Row-Major.
  gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
  gsl_matrix_transpose_memcpy(Tt_orig, T);

  // Pre-generate all permutation indices
  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  for(size_t k = 0; k < n; k++) temp_idx[k] = (int)k;

  srand(0);
  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, (int)n);
    memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
  }
  free(temp_idx);

  R_xlen_t total_elements = (R_xlen_t)p * (R_xlen_t)m;

  // Initialize accumulators
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i = 0; i < total_elements; i++) {
    pvalue_vec[i] = 0.0;
    se_vec[i] = 0.0;
    zscore_vec[i] = 0.0;
  }

  int actual_threads = 1;
  #ifdef _OPENMP
    #pragma omp parallel
    {
      #pragma omp single
      actual_threads = omp_get_num_threads();
    }
  #endif

  // Per-thread reduction buffers
  double *thread_sum = (double*)calloc((size_t)actual_threads * total_elements, sizeof(double));
  double *thread_sum_sq = (double*)calloc((size_t)actual_threads * total_elements, sizeof(double));
  double *thread_count = (double*)calloc((size_t)actual_threads * total_elements, sizeof(double));

  // Parallel permutation loop
  #pragma omp parallel
  {
    int tid = 0;
    #ifdef _OPENMP
      tid = omp_get_thread_num();
    #endif

    // Pointers to thread-local accumulators
    double *my_sum = &thread_sum[(size_t)tid * total_elements];
    double *my_sum_sq = &thread_sum_sq[(size_t)tid * total_elements];
    double *my_count = &thread_count[(size_t)tid * total_elements];

    // Thread-local working matrices
    gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

    // Helpers for memcpy
    double *src_base = Tt_orig->data;
    double *dst_base = Tt_perm->data;
    size_t row_size_bytes = p * sizeof(double);

    #pragma omp for schedule(dynamic, 16)
    for(int i_rand = 0; i_rand < nrand; i_rand++) 
    {
      int *p_idx = &perm_table[(size_t)i_rand * n];
      
      // OPTIMIZED PERMUTATION: Row-wise memcpy
      // Tt_perm[row i] = Tt_orig[row p_idx[i]]
      for(size_t i = 0; i < n; i++) {
        int src_row_idx = p_idx[i];
        memcpy(dst_base + (i * p), 
               src_base + (src_row_idx * p), 
               row_size_bytes);
      }

      // beta_rand = (Tt_perm)^T @ (Yt)^T
      gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Tt_perm, Yt, 0.0, beta_rand);

      // Accumulate
      for(size_t r = 0; r < p; r++) {
        for(size_t c = 0; c < m; c++) {
          R_xlen_t idx = (R_xlen_t)r * m + c;
          double b_rnd = gsl_matrix_get(beta_rand, r, c);
          double b_obs = beta_vec[idx];
          
          my_sum[idx] += b_rnd;
          my_sum_sq[idx] += b_rnd * b_rnd;
          if(fabs(b_rnd) >= fabs(b_obs)) {
            my_count[idx] += 1.0;
          }
        }
      }
    }

    gsl_matrix_free(Tt_perm);
    gsl_matrix_free(beta_rand);
  }

  // Reduce results
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i = 0; i < total_elements; i++) {
    double sum_val = 0.0, sum_sq_val = 0.0, count_val = 0.0;
    for(int t = 0; t < actual_threads; t++) {
      size_t offset = (size_t)t * total_elements + i;
      sum_val += thread_sum[offset];
      sum_sq_val += thread_sum_sq[offset];
      count_val += thread_count[offset];
    }
    zscore_vec[i] = sum_val;
    se_vec[i] = sum_sq_val;
    pvalue_vec[i] = count_val;
  }

  free(thread_sum);
  free(thread_sum_sq);
  free(thread_count);
  free(perm_table);

  // Finalize statistics
  double inv_nrand = 1.0 / nrand;
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i = 0; i < total_elements; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean_rand = zscore_vec[i] * inv_nrand;
    double var_rand = (se_vec[i] * inv_nrand) - (mean_rand * mean_rand);
    double se_rand = sqrt(var_rand > 0 ? var_rand : 0);
    se_vec[i] = se_rand;
    zscore_vec[i] = (se_rand > 1e-12) ? (beta_vec[i] - mean_rand) / se_rand : 0.0;
  }

  gsl_matrix_free(Tt_orig);
  gsl_matrix_free(I_mat);
  gsl_matrix_free(T);
  gsl_matrix_partial_free(Xt);
  gsl_matrix_partial_free(Yt);
  gsl_matrix_partial_free(beta);
}

SEXP ridgeRegTperm_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp, SEXP ncores_sexp) 
{
  SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
  size_t n = (size_t)INTEGER(x_dim)[0];  // genes
  size_t p = (size_t)INTEGER(x_dim)[1];  // proteins

  SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
  size_t m = (size_t)INTEGER(y_dim)[1];  // samples
  
  R_xlen_t total_len = (R_xlen_t)p * (R_xlen_t)m;

  SEXP beta_sexp   = PROTECT(allocVector(REALSXP, total_len));
  SEXP se_sexp     = PROTECT(allocVector(REALSXP, total_len));
  SEXP zscore_sexp = PROTECT(allocVector(REALSXP, total_len));
  SEXP pvalue_sexp = PROTECT(allocVector(REALSXP, total_len));

  ridgeRegTperm_core(
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

/* =========================================================
   REGISTRATION
   ========================================================= */

static const R_CMethodDef cMethods[] = {
  {"ridgeReg", (DL_FUNC) &ridgeReg, 11},
  {"ridgeRegTperm_single", (DL_FUNC) &ridgeRegTperm_single, 11},
  {NULL, NULL, 0}
};

static const R_CallMethodDef callMethods[] = {
  {"ridgeRegFast_interface",  (DL_FUNC) &ridgeRegFast_interface,  5},
  {"ridgeRegTperm_interface", (DL_FUNC) &ridgeRegTperm_interface, 5},
  {NULL, NULL, 0}
};

void R_init_RidgeR(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
