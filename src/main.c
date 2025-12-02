#include <R.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <string.h> // for memcpy
#include <omp.h>    // for OpenMP
#include <time.h>   // for seed
#include <math.h>   // for sqrt, fabs

// Tuning parameter: Number of permutations to process in one matrix multiplication
// 128 is usually a sweet spot for cache efficiency.
#define BLOCK_SIZE 128 

/* =========================================================
   HELPER FUNCTIONS
   ========================================================= */

// Wrapper helper: Maps R vectors to GSL matrix structs
gsl_matrix *RVectorObject_to_gsl_matrix(double *vec, size_t nr, size_t nc)
{
  gsl_block *b;
  gsl_matrix *r;

  b = (gsl_block*)malloc(sizeof(gsl_block));
  r = (gsl_matrix*)malloc(sizeof(gsl_matrix));

  r->size1 = nr;
  r->tda = r->size2 = nc;
  r->owner = 1; // Mark as owner of the struct (but not the data)

  b->data = r->data = vec;
  r->block = b;
  b->size = r->size1 * r->size2;

  return r;
}

// Custom free helper for the wrapper
void gsl_matrix_partial_free(gsl_matrix *x)
{
  if(x) {
    if(x->block) free(x->block);
    free(x);
  }
}

// Legacy Shuffle: Used by Old Version AND Fast Version (for pre-calculation)
// to ensure results are identical.
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
   OLD VERSION (Legacy, Slow, Single Thread)
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

  // I = (X'X + lambda*I)
  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);

  // I = (X'X + lambda*I)^-1
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);

  // T = (X'X + lambda*I)^-1 * X'
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);

  // beta = T * Y
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

  // Permutation
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
   NEW VERSION (Fast, OpenMP, Batched + Deterministic RNG)
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

  // 3. Pre-calculate Permutations (Deterministic & Identical to Old)
  // We allocate a table to store all permutation indices. 
  // We fill this sequentially using the EXACT logic as the Old function.
  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  
  // Initialize to 0..n-1
  int k;
  for(k=0; k<n; k++) temp_idx[k] = k;

  // Reset Seed to match the Old function
  srand(0);

  // Generate shuffled indices
  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, n);
    // Store in the table
    memcpy(&perm_table[i_rand * n], temp_idx, n * sizeof(int));
  }
  free(temp_idx);

  // 4. Batched Parallel Execution
  gsl_matrix_set_zero(pvalue);
  gsl_matrix_set_zero(aver_sq);

  #pragma omp parallel
  {
    // Thread-local accumulators
    gsl_matrix *pvalue_priv = gsl_matrix_calloc(p, m);
    gsl_matrix *aver_priv = gsl_matrix_calloc(p, m);
    gsl_matrix *aver_sq_priv = gsl_matrix_calloc(p, m);
    
    // Batch Buffers:
    // Y_block stacks 'BLOCK_SIZE' number of Y matrices horizontally.
    gsl_matrix *Y_block = gsl_matrix_alloc(n, m * BLOCK_SIZE);
    
    // Beta_block stores the result of T * Y_block
    gsl_matrix *Beta_block = gsl_matrix_alloc(p, m * BLOCK_SIZE);
    
    int b_start;
    
    // Iterate through batches in parallel
    #pragma omp for schedule(dynamic)
    for(b_start = 0; b_start < nrand; b_start += BLOCK_SIZE)
    {
      // Calculate how many permutations are in this batch (handle end of loop)
      int current_batch_sz = (b_start + BLOCK_SIZE > nrand) ? (nrand - b_start) : BLOCK_SIZE;
      
      // Step A: Construct Y_block using Pre-calculated Permutations
      int i_perm, row;
      for(i_perm = 0; i_perm < current_batch_sz; i_perm++) {
        
        // Find the specific permutation indices for this iteration
        int global_perm_id = b_start + i_perm;
        int *current_perm_idx = &perm_table[global_perm_id * n];

        // Copy permuted Y rows into the block matrix
        for(row = 0; row < n; row++) {
           memcpy(Y_block->data + (row * Y_block->tda) + (i_perm * m),
                  Y->data + (current_perm_idx[row] * Y->tda),
                  m * sizeof(double));
        }
      }

      // Step B: Batched Matrix Multiplication (High Performance)
      // View into the valid part of the buffers
      gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block, 0, 0, n, m * current_batch_sz);
      gsl_matrix_view B_sub = gsl_matrix_submatrix(Beta_block, 0, 0, p, m * current_batch_sz);
      
      // Beta_block = T * Y_block
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);

      // Step C: Accumulate Statistics
      int r, c;
      for(i_perm = 0; i_perm < current_batch_sz; i_perm++) {
        for(r = 0; r < p; r++) {
          for(c = 0; c < m; c++) {
             // Access randomized beta
             double b_rnd = gsl_matrix_get(&B_sub.matrix, r, (i_perm * m) + c);
             // Access observed beta
             double b_obs = beta->data[r * m + c];

             if(fabs(b_rnd) >= fabs(b_obs)) {
                pvalue_priv->data[r * m + c] += 1.0;
             }
             aver_priv->data[r * m + c] += b_rnd;
             aver_sq_priv->data[r * m + c] += (b_rnd * b_rnd);
          }
        }
      }
    }

    // Reduce thread-local results to global matrices
    #pragma omp critical
    {
      gsl_matrix_add(pvalue, pvalue_priv);
      gsl_matrix_add(aver, aver_priv);
      gsl_matrix_add(aver_sq, aver_sq_priv);
    }

    // Free thread-local memory
    gsl_matrix_free(Y_block);
    gsl_matrix_free(Beta_block);
    gsl_matrix_free(pvalue_priv);
    gsl_matrix_free(aver_priv);
    gsl_matrix_free(aver_sq_priv);
  }

  // Done with permutation table
  free(perm_table);

  // 5. Final Statistics Calculation (Matches Old Code Logic)
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
