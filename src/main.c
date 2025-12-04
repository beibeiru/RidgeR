#include <R.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <string.h> // for memcpy
#include <time.h>   // for seed
#include <math.h>   // for sqrt, fabs

#ifdef _OPENMP
  #include <omp.h>  // for OpenMP
#endif

// --- TUNING ---
// Number of samples (columns) a thread processes at one time.
// 64 is small enough to fit in cache, large enough for BLAS efficiency.
// This keeps memory usage CONSTANT per thread, regardless of m.
#define SAMPLE_STRIP_SIZE 64 

// Batch of permutations. 
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
   NEW VERSION (Safe for Huge Data)
   ========================================================= */
void ridgeRegFast(
  double *X_vec, double *Y_vec,
  int *n_pt, int *p_pt, int *m_pt,
  double *lambda_pt, double *nrand_pt,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec,
  int *nthreads_pt
)
{
  // Use size_t for all dimensions to prevent integer overflow
  size_t n = (size_t)*n_pt;
  size_t p = (size_t)*p_pt;
  size_t m = (size_t)*m_pt;
  int nrand = (int)*nrand_pt;
  double lambda = *lambda_pt;
  int num_threads = *nthreads_pt;

  if (num_threads > 0) {
       #ifdef _OPENMP
           omp_set_num_threads(num_threads);
       #endif
  }

  gsl_matrix *X = RVectorObject_to_gsl_matrix(X_vec, n, p);
  gsl_matrix *Y = RVectorObject_to_gsl_matrix(Y_vec, n, m);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);
  
  // NOTE: For huge m, we do NOT allocate full accumulators in GSL
  // We use the raw double pointers passed from R directly (se_vec, etc.)
  
  gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);
  
  // 1. Projection Matrix
  gsl_matrix_set_identity(I_mat);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I_mat);
  gsl_linalg_cholesky_decomp(I_mat);
  gsl_linalg_cholesky_invert(I_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I_mat, X, 0.0, T);

  // 2. Real Beta
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

  // 3. Pre-calculate Permutations
  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  
  size_t k;
  for(k=0; k<n; k++) temp_idx[k] = (int)k;

  srand(0);

  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, (int)n);
    memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
  }
  free(temp_idx);

  // Initialize outputs to 0
  // Note: Using parallel loop for initialization as size could be huge (1M * 1000)
  size_t total_elements = p * m;
  #pragma omp parallel for
  for(size_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = 0.0;
    se_vec[i] = 0.0;     // aver_sq
    zscore_vec[i] = 0.0; // aver
  }
  // Reuse pointers for clarity: zscore_vec is used for 'aver' temporarily

  // 4. Parallelize over SAMPLES (Columns)
  // Each thread takes a disjoint set of columns. No memory conflicts.
  #pragma omp parallel
  {
    // Small buffers, fixed size. Safe for any m.
    // Dimensions: n x (SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE)
    gsl_matrix *Y_block = gsl_matrix_alloc(n, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    gsl_matrix *Beta_block = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    
    // Iterate over chunks of columns (Samples)
    size_t col_start;
    #pragma omp for schedule(dynamic)
    for(col_start = 0; col_start < m; col_start += SAMPLE_STRIP_SIZE) 
    {
      size_t current_cols = (col_start + SAMPLE_STRIP_SIZE > m) ? (m - col_start) : SAMPLE_STRIP_SIZE;
      
      // Iterate over batches of permutations
      int b_start;
      for(b_start = 0; b_start < nrand; b_start += PERM_BATCH_SIZE)
      {
        int current_perms = (b_start + PERM_BATCH_SIZE > nrand) ? (nrand - b_start) : PERM_BATCH_SIZE;

        // A. Construct Y_block (Stacking permutations side-by-side)
        int i_perm;
        size_t row;
        
        for(i_perm = 0; i_perm < current_perms; i_perm++) {
           int *current_perm_idx = &perm_table[(size_t)(b_start + i_perm) * n];
           
           for(row = 0; row < n; row++) {
             // Copy logic: Get permuted ROW from Y, copy 'current_cols' elements
             // Y->data is row-major? No, R matrices are usually COL-major, 
             // but our wrapper says size1=n, size2=m. 
             // IMPORTANT: R stores matrices as Column-Major.
             // But we passed t(X) and t(Y) from R, so in C they look Row-Major.
             // So Y->data is linear. row 'r' starts at Y->data + r*m.
             
             // Offset in Y_block: (row * tda) + (i_perm * current_cols)
             // Offset in Y:       (permuted_row * tda) + col_start
             
             memcpy(Y_block->data + (row * Y_block->tda) + (i_perm * current_cols),
                    Y->data + ((size_t)current_perm_idx[row] * Y->tda) + col_start,
                    current_cols * sizeof(double));
           }
        }

        // B. Matrix Multiply
        gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block, 0, 0, n, current_perms * current_cols);
        gsl_matrix_view B_sub = gsl_matrix_submatrix(Beta_block, 0, 0, p, current_perms * current_cols);
        
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);

        // C. Update Stats (Direct write to global arrays, no contention)
        size_t r, c_local;
        for(i_perm = 0; i_perm < current_perms; i_perm++) {
          for(r = 0; r < p; r++) {
            for(c_local = 0; c_local < current_cols; c_local++) {
              
              size_t global_col_idx = col_start + c_local;
              size_t global_idx = r * m + global_col_idx; // Linear index in output
              
              double b_rnd = gsl_matrix_get(&B_sub.matrix, r, (i_perm * current_cols) + c_local);
              double b_obs = beta_vec[global_idx]; 

              if(fabs(b_rnd) >= fabs(b_obs)) {
                pvalue_vec[global_idx] += 1.0;
              }
              zscore_vec[global_idx] += b_rnd;           // Sum
              se_vec[global_idx]     += (b_rnd * b_rnd); // Sum Sq
            }
          }
        }
      }
    }

    gsl_matrix_free(Y_block);
    gsl_matrix_free(Beta_block);
  }

  free(perm_table);

  // 5. Finalize (Parallelized)
  double inv_nrand = 1.0 / nrand;
  
  #pragma omp parallel for
  for(size_t i=0; i<total_elements; i++) {
    // P-value
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);

    // Mean and SE
    double mean_rand = zscore_vec[i] * inv_nrand;
    double mean_sq_rand = se_vec[i] * inv_nrand;
    
    double var_rand = mean_sq_rand - (mean_rand * mean_rand);
    if(var_rand < 0) var_rand = 0;
    double se_rand = sqrt(var_rand);
    
    // Store back
    se_vec[i] = se_rand; // This is now SE, not SumSq
    
    // Z-score
    double obs = beta_vec[i];
    if(se_rand > 0) {
      zscore_vec[i] = (obs - mean_rand) / se_rand;
    } else {
      zscore_vec[i] = 0.0;
    }
  }

  gsl_matrix_free(I_mat);
  gsl_matrix_free(T);
  // Note: 'aver' was not allocated, we used zscore_vec
  // gsl_matrix_free(aver); 

  gsl_matrix_partial_free(X);
  gsl_matrix_partial_free(Y);
  gsl_matrix_partial_free(beta);
  // gsl_matrix_partial_free for se/zscore/pvalue is not needed as we didn't use structs for accumulators
}
