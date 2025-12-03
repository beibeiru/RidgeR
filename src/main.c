#include <R.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <string.h> // for memcpy
#include <omp.h>    // for OpenMP
#include <time.h>   // for seed
#include <math.h>   // for sqrt, fabs

// --- MEMORY SAFETY LIMITS ---
// We limit the buffer width to ensure we never exhaust RAM.
// 2048 columns * 20,000 rows * 8 bytes ~= 320 MB per thread.
// This is safe even for 32+ cores on standard servers.
#define MAX_BUFFER_COLS 2048 

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
   OLD VERSION
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

  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

  srand(0);
  for(i=0;i<n;i++) array_index[i] = i;

  gsl_matrix_set_zero(aver);
  gsl_matrix_set_zero(aver_sq);
  gsl_matrix_set_zero(pvalue);

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
   NEW VERSION (Double Tiling: Robust for m=1 to m=1,000,000)
   ========================================================= */
void ridgeRegFast(
  double *X_vec, double *Y_vec,
  int *n_pt, int *p_pt, int *m_pt,
  double *lambda_pt, double *nrand_pt,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int *nthreads_pt
)
{
  int n = *n_pt, p = *p_pt, m = *m_pt, nrand = (int)*nrand_pt;
  double lambda = *lambda_pt;
  int num_threads = *nthreads_pt;

  if (num_threads > 0) omp_set_num_threads(num_threads);

  gsl_matrix *X = RVectorObject_to_gsl_matrix(X_vec, n, p);
  gsl_matrix *Y = RVectorObject_to_gsl_matrix(Y_vec, n, m);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);
  gsl_matrix *aver_sq = RVectorObject_to_gsl_matrix(se_vec, p, m);
  gsl_matrix *zscore = RVectorObject_to_gsl_matrix(zscore_vec, p, m);
  gsl_matrix *pvalue = RVectorObject_to_gsl_matrix(pvalue_vec, p, m);

  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);
  gsl_matrix *aver = gsl_matrix_calloc(p, m);

  // 1. Projection Matrix
  gsl_matrix_set_identity(I_mat);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I_mat);
  gsl_linalg_cholesky_decomp(I_mat);
  gsl_linalg_cholesky_invert(I_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I_mat, X, 0.0, T);

  // 2. Real Beta
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

  // 3. Pre-calculate Permutations (Deterministic)
  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  
  int k;
  for(k=0; k<n; k++) temp_idx[k] = k;

  srand(0);

  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, n);
    memcpy(&perm_table[i_rand * n], temp_idx, n * sizeof(int));
  }
  free(temp_idx);

  // 4. Batched Parallel Execution
  // We determine Strip Size and Batch Size to ensure we never exceed MAX_BUFFER_COLS
  
  // Decide how many columns of Y to process at once (Strip Size)
  // If m is huge, we restrict this to allow room for permutation stacking.
  // We aim to process at least 8 permutations at once for BLAS efficiency.
  int strip_size = m;
  int min_perms = 8;
  
  if (strip_size * min_perms > MAX_BUFFER_COLS) {
    strip_size = MAX_BUFFER_COLS / min_perms;
  }
  if (strip_size < 1) strip_size = 1;

  // Now determine how many permutations fit in the remaining space
  int batch_size = MAX_BUFFER_COLS / strip_size;
  if (batch_size < 1) batch_size = 1;

  gsl_matrix_set_zero(pvalue);
  gsl_matrix_set_zero(aver_sq);

  #pragma omp parallel
  {
    // Thread-local output accumulators (Full Size)
    gsl_matrix *pvalue_priv = gsl_matrix_calloc(p, m);
    gsl_matrix *aver_priv = gsl_matrix_calloc(p, m);
    gsl_matrix *aver_sq_priv = gsl_matrix_calloc(p, m);
    
    // Safety Buffers (Fixed Max Size)
    // We allocate these based on MAX_BUFFER_COLS, not m.
    gsl_matrix *Y_block = gsl_matrix_alloc(n, strip_size * batch_size);
    gsl_matrix *Beta_block = gsl_matrix_alloc(p, strip_size * batch_size);
    
    // OUTER LOOP: Tile over Samples (Columns of Y)
    // This handles huge m by breaking it into chunks.
    int c_start;
    #pragma omp for schedule(dynamic)
    for(c_start = 0; c_start < m; c_start += strip_size)
    {
       int current_strip_w = (c_start + strip_size > m) ? (m - c_start) : strip_size;
       
       // INNER LOOP: Tile over Permutations
       // This handles nrand using efficient batches
       int b_start;
       for(b_start = 0; b_start < nrand; b_start += batch_size) 
       {
           int current_batch_sz = (b_start + batch_size > nrand) ? (nrand - b_start) : batch_size;
           
           // Step A: Fill Y_block
           // We stack 'current_batch_sz' copies of the current Y strip side-by-side
           int i_perm, row, col;
           
           for(i_perm = 0; i_perm < current_batch_sz; i_perm++) {
             int global_perm_id = b_start + i_perm;
             int *current_perm_idx = &perm_table[global_perm_id * n];
             
             // Copy logic
             // Source: Y (starts at column c_start, width current_strip_w)
             // Dest: Y_block (starts at column i_perm * current_strip_w)
             for(row = 0; row < n; row++) {
                // We copy a row segment (width = current_strip_w) via memcpy
                // But we must account for shuffling. 
                // The source ROW is permuted. The source COLUMNS are c_start to c_start + w.
                memcpy(Y_block->data + (row * Y_block->tda) + (i_perm * current_strip_w),
                       Y->data + (current_perm_idx[row] * Y->tda) + c_start,
                       current_strip_w * sizeof(double));
             }
           }
           
           // Step B: Batched Matrix Multiply
           int total_block_width = current_batch_sz * current_strip_w;
           gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block, 0, 0, n, total_block_width);
           gsl_matrix_view B_sub = gsl_matrix_submatrix(Beta_block, 0, 0, p, total_block_width);
           
           gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);
           
           // Step C: Accumulate
           int r, c_local;
           for(i_perm = 0; i_perm < current_batch_sz; i_perm++) {
             for(r = 0; r < p; r++) {
               for(c_local = 0; c_local < current_strip_w; c_local++) {
                  
                  // Map local coordinates back to Global Y coordinates
                  int global_c = c_start + c_local;
                  
                  double b_rnd = gsl_matrix_get(&B_sub.matrix, r, (i_perm * current_strip_w) + c_local);
                  double b_obs = beta->data[r * m + global_c];
                  
                  // Update stats at global coordinate
                  if(fabs(b_rnd) >= fabs(b_obs)) {
                     pvalue_priv->data[r * m + global_c] += 1.0;
                  }
                  aver_priv->data[r * m + global_c] += b_rnd;
                  aver_sq_priv->data[r * m + global_c] += (b_rnd * b_rnd);
               }
             }
           }
       }
    }

    #pragma omp critical
    {
      gsl_matrix_add(pvalue, pvalue_priv);
      gsl_matrix_add(aver, aver_priv);
      gsl_matrix_add(aver_sq, aver_sq_priv);
    }

    gsl_matrix_free(Y_block);
    gsl_matrix_free(Beta_block);
    gsl_matrix_free(pvalue_priv);
    gsl_matrix_free(aver_priv);
    gsl_matrix_free(aver_sq_priv);
  }

  free(perm_table);

  // 5. Final Stats
  double inv_nrand = 1.0 / nrand;

  gsl_matrix_add_constant(pvalue, 1.0);
  gsl_matrix_scale(pvalue, 1.0 / (nrand + 1));

  gsl_matrix_scale(aver, inv_nrand);
  gsl_matrix_scale(aver_sq, inv_nrand);

  gsl_matrix_memcpy(zscore, beta);
  gsl_matrix_sub(zscore, aver);

  gsl_matrix_mul_elements(aver, aver);
  gsl_matrix_sub(aver_sq, aver);

  for(int x = 0; x < p * m; x++) {
    if(aver_sq->data[x] < 0) aver_sq->data[x] = 0;
    aver_sq->data[x] = sqrt(aver_sq->data[x]);
  }

  gsl_matrix_div_elements(zscore, aver_sq);

  gsl_matrix_free(I_mat);
  gsl_matrix_free(T);
  gsl_matrix_free(aver);
  gsl_matrix_partial_free(X);
  gsl_matrix_partial_free(Y);
  gsl_matrix_partial_free(beta);
  gsl_matrix_partial_free(aver_sq);
  gsl_matrix_partial_free(zscore);
  gsl_matrix_partial_free(pvalue);
}
