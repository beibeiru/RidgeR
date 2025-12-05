#include <R.h>
#include <Rinternals.h> // Required for .Call
#include <R_ext/Rdynload.h> // For registration
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
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
  r->owner = 1; // We own the struct, but NOT the data
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

// Internal logic for the new function
void ridgeRegFast_core(
  double *X_vec, double *Y_vec,
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

  // Use R_xlen_t for total elements to support > 2 billion
  R_xlen_t total_elements = (R_xlen_t)p * (R_xlen_t)m;

  gsl_matrix *X = RVectorObject_to_gsl_matrix(X_vec, n, p);
  gsl_matrix *Y = RVectorObject_to_gsl_matrix(Y_vec, n, m);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);
  
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

  // 3. Permutations Setup
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
  // Using R_xlen_t for loop variable to handle long vectors
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = 0.0;
    se_vec[i] = 0.0;     
    zscore_vec[i] = 0.0; 
  }

  // 4. Parallel Processing
  #pragma omp parallel
  {
    gsl_matrix *Y_block = gsl_matrix_alloc(n, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    gsl_matrix *Beta_block = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    
    #pragma omp for schedule(dynamic)
    for(size_t col_start = 0; col_start < m; col_start += SAMPLE_STRIP_SIZE) 
    {
      size_t current_cols = (col_start + SAMPLE_STRIP_SIZE > m) ? (m - col_start) : SAMPLE_STRIP_SIZE;
      
      for(int b_start = 0; b_start < nrand; b_start += PERM_BATCH_SIZE)
      {
        int current_perms = (b_start + PERM_BATCH_SIZE > nrand) ? (nrand - b_start) : PERM_BATCH_SIZE;

        // Construct Permuted Y Block
        for(int i_perm = 0; i_perm < current_perms; i_perm++) {
           int *current_perm_idx = &perm_table[(size_t)(b_start + i_perm) * n];
           for(size_t row = 0; row < n; row++) {
             memcpy(Y_block->data + (row * Y_block->tda) + ((size_t)i_perm * current_cols),
                    Y->data + ((size_t)current_perm_idx[row] * Y->tda) + col_start,
                    current_cols * sizeof(double));
           }
        }

        gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block, 0, 0, n, current_perms * current_cols);
        gsl_matrix_view B_sub = gsl_matrix_submatrix(Beta_block, 0, 0, p, current_perms * current_cols);
        
        // Matrix Multiply
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);

        // Aggregate Stats
        for(int i_perm = 0; i_perm < current_perms; i_perm++) {
          for(size_t r = 0; r < p; r++) {
            for(size_t c_local = 0; c_local < current_cols; c_local++) {
              
              R_xlen_t global_idx = (R_xlen_t)r * m + (col_start + c_local);
              
              double b_rnd = gsl_matrix_get(&B_sub.matrix, r, ((size_t)i_perm * current_cols) + c_local);
              double b_obs = beta_vec[global_idx]; 

              if(fabs(b_rnd) >= fabs(b_obs)) {
                pvalue_vec[global_idx] += 1.0;
              }
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

  double inv_nrand = 1.0 / nrand;
  #pragma omp parallel for schedule(static)
  for(R_xlen_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean_rand = zscore_vec[i] * inv_nrand;
    double mean_sq_rand = se_vec[i] * inv_nrand;
    double var_rand = mean_sq_rand - (mean_rand * mean_rand);
    if(var_rand < 0) var_rand = 0;
    double se_rand = sqrt(var_rand);
    se_vec[i] = se_rand;
    double obs = beta_vec[i];
    zscore_vec[i] = (se_rand > 1e-12) ? (obs - mean_rand) / se_rand : 0.0;
  }

  gsl_matrix_free(I_mat); gsl_matrix_free(T);
  gsl_matrix_partial_free(X); gsl_matrix_partial_free(Y); gsl_matrix_partial_free(beta);
}

// .Call Interface (Required for Long Vectors)
SEXP ridgeRegFast_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp, SEXP ncores_sexp) 
{
    double lambda = asReal(lambda_sexp);
    int nrand = asInteger(nrand_sexp);
    int ncores = asInteger(ncores_sexp);

    // Get dims. Note: Inputs are Transposed matrices from R.
    // X is (P x N), Y is (M x N)
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t p = (size_t)INTEGER(x_dim)[0]; 
    size_t n = (size_t)INTEGER(x_dim)[1];

    SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
    size_t m = (size_t)INTEGER(y_dim)[0];
    
    // Allocate Huge Vectors (Long Vector support)
    R_xlen_t total_len = (R_xlen_t)p * (R_xlen_t)m;

    SEXP beta_sexp   = PROTECT(allocVector(REALSXP, total_len));
    SEXP se_sexp     = PROTECT(allocVector(REALSXP, total_len));
    SEXP zscore_sexp = PROTECT(allocVector(REALSXP, total_len));
    SEXP pvalue_sexp = PROTECT(allocVector(REALSXP, total_len));

    ridgeRegFast_core(
        REAL(X_sexp), REAL(Y_sexp), n, p, m, lambda, nrand,
        REAL(beta_sexp), REAL(se_sexp), REAL(zscore_sexp), REAL(pvalue_sexp),
        ncores
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


