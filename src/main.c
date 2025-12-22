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
  int i, j, t;
  for (i = 0; i < n - 1; i++) {
    j = i + rand() / (RAND_MAX / (n - i) + 1);
    t = array[j];
    array[j] = array[i];
    array[i] = t;
  }
}

/* =========================================================
   VERSION 1: OLD (Single-threaded, Y row permutation)
   ========================================================= */
void ridgeReg_old_core(
  double *X_ptr, double *Y_ptr, size_t n, size_t p, size_t m,
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec
)
{
  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n);
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);
  gsl_matrix *Y_rand_t = gsl_matrix_alloc(m, n);
  gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
  
  double *sum = (double*)calloc(p * m, sizeof(double));
  double *sum_sq = (double*)calloc(p * m, sizeof(double));

  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I);
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, Xt, 0.0, T);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  int *array_index = (int*)malloc(n * sizeof(int));
  for(size_t i=0; i<n; i++) array_index[i] = (int)i;

  srand(0);
  for(int i_rand=0; i_rand<nrand; i_rand++)
  {
    shuffle(array_index, (int)n);
    for(size_t j=0; j<n; j++){
      gsl_vector_const_view t_col = gsl_matrix_const_column(Yt, (size_t)array_index[j]);
      // FIXED: Use gsl_matrix_set_col instead of gsl_matrix_set_column
      gsl_matrix_set_col(Y_rand_t, j, &t_col.vector);
    }
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand_t, 0.0, beta_rand);

    for(size_t k=0; k < p*m; k++) {
      double br = beta_rand->data[k];
      if(fabs(br) >= fabs(beta_vec[k])) pvalue_vec[k]++;
      sum[k] += br;
      sum_sq[k] += br * br;
    }
  }

  double inv_nr = 1.0 / nrand;
  for(size_t i=0; i < p*m; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean = sum[i] * inv_nr;
    double var = (sum_sq[i] * inv_nr) - (mean * mean);
    double sd = sqrt(var > 0 ? var : 0);
    se_vec[i] = sd;
    zscore_vec[i] = (sd > 1e-12) ? (beta_vec[i] - mean) / sd : 0.0;
  }

  gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_rand_t);
  gsl_matrix_free(beta_rand); free(array_index); free(sum); free(sum_sq);
  gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta);
}

/* =========================================================
   VERSION 1.5: OLD2 (Single-threaded, T permutation)
   ========================================================= */
void ridgeRegTperm_old_core(
  double *X_ptr, double *Y_ptr, size_t n, size_t p, size_t m,
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec
)
{
  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n);
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);
  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I);
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, Xt, 0.0, T);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
  gsl_matrix_transpose_memcpy(Tt_orig, T);
  gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);
  gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

  int *temp_idx = (int*)malloc(n * sizeof(int));
  for(size_t k=0; k<n; k++) temp_idx[k] = (int)k;
  
  double *sum = (double*)calloc(p * m, sizeof(double));
  double *sum_sq = (double*)calloc(p * m, sizeof(double));

  srand(0);
  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, (int)n);
    for(size_t i=0; i<n; i++) {
       memcpy(Tt_perm->data + ((size_t)temp_idx[i] * p), Tt_orig->data + (i * p), p * sizeof(double));
    }
    gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Tt_perm, Yt, 0.0, beta_rand);
    for(size_t k=0; k < p*m; k++) {
      double br = beta_rand->data[k];
      if(fabs(br) >= fabs(beta_vec[k])) pvalue_vec[k]++;
      sum[k] += br; sum_sq[k] += br * br;
    }
  }

  double inv_nr = 1.0 / nrand;
  for(size_t i=0; i < p*m; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean = sum[i] * inv_nr;
    double sd = sqrt(fmax(0, (sum_sq[i] * inv_nr) - (mean * mean)));
    se_vec[i] = sd;
    zscore_vec[i] = (sd > 1e-12) ? (beta_vec[i] - mean) / sd : 0.0;
  }

  gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Tt_orig); 
  gsl_matrix_free(Tt_perm); gsl_matrix_free(beta_rand);
  free(temp_idx); free(sum); free(sum_sq);
  gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta);
}

/* =========================================================
   VERSION 2: NEW (Multi-threaded, Y row permutation)
   ========================================================= */
void ridgeRegFast_core(
  double *X_ptr, double *Y_ptr, size_t n, size_t p, size_t m,
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec, int num_threads
)
{
  #ifdef _OPENMP
    if (num_threads > 0) omp_set_num_threads(num_threads);
  #endif

  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n);
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);
  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I);
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, Xt, 0.0, T);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  for(size_t k=0; k<n; k++) temp_idx[k] = (int)k;
  srand(0);
  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, (int)n);
    memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
  }
  free(temp_idx);

  size_t total_elements = p * m;
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = 0.0; se_vec[i] = 0.0; zscore_vec[i] = 0.0; 
  }

  #pragma omp parallel
  {
    gsl_matrix *Y_block = gsl_matrix_alloc(n, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    gsl_matrix *B_block = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
    
    #pragma omp for schedule(dynamic)
    for(size_t samp_start = 0; samp_start < m; samp_start += SAMPLE_STRIP_SIZE) 
    {
      size_t cur_s = (samp_start + SAMPLE_STRIP_SIZE > m) ? (m - samp_start) : SAMPLE_STRIP_SIZE;
      for(int b_start = 0; b_start < nrand; b_start += PERM_BATCH_SIZE)
      {
        int cur_p = (b_start + PERM_BATCH_SIZE > nrand) ? (nrand - b_start) : PERM_BATCH_SIZE;
        for(int ip = 0; ip < cur_p; ip++) {
          int *p_idx = &perm_table[(size_t)(b_start + ip) * n];
          for(size_t sl = 0; sl < cur_s; sl++) {
            double *src_row = gsl_matrix_ptr(Yt, samp_start + sl, 0);
            size_t dest_col = (ip * cur_s) + sl;
            for(size_t g = 0; g < n; g++) gsl_matrix_set(Y_block, g, dest_col, src_row[p_idx[g]]);
          }
        }
        gsl_matrix_view Ys = gsl_matrix_submatrix(Y_block, 0, 0, n, cur_p * cur_s);
        gsl_matrix_view Bs = gsl_matrix_submatrix(B_block, 0, 0, p, cur_p * cur_s);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Ys.matrix, 0.0, &Bs.matrix);

        for(int ip = 0; ip < cur_p; ip++) {
          for(size_t r = 0; r < p; r++) {
            for(size_t sl = 0; sl < cur_s; sl++) {
              size_t idx = r * m + (samp_start + sl);
              double br = gsl_matrix_get(&Bs.matrix, r, (ip * cur_s) + sl);
              if(fabs(br) >= fabs(beta_vec[idx])) pvalue_vec[idx] += 1.0;
              zscore_vec[idx] += br; se_vec[idx] += (br * br);
            }
          }
        }
      }
    }
    gsl_matrix_free(Y_block); gsl_matrix_free(B_block);
  }

  double inv_nr = 1.0 / nrand;
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<total_elements; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean = zscore_vec[i] * inv_nr;
    double sd = sqrt(fmax(0, (se_vec[i] * inv_nr) - (mean * mean)));
    se_vec[i] = sd;
    zscore_vec[i] = (sd > 1e-12) ? (beta_vec[i] - mean) / sd : 0.0;
  }
  free(perm_table); gsl_matrix_free(I); gsl_matrix_free(T);
  gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta);
}

/* =========================================================
   VERSION 3: NEW2 (Multi-threaded, T permutation)
   ========================================================= */
void ridgeRegTperm_core(
  double *X_ptr, double *Y_ptr, size_t n, size_t p, size_t m,
  double lambda, int nrand,
  double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec, int num_threads
)
{
  #ifdef _OPENMP
    if (num_threads > 0) omp_set_num_threads(num_threads);
  #endif

  gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n);
  gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n);
  gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  gsl_matrix *I = gsl_matrix_alloc(p, p);
  gsl_matrix *T = gsl_matrix_alloc(p, n);
  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I);
  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, Xt, 0.0, T);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta);

  gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
  gsl_matrix_transpose_memcpy(Tt_orig, T);

  int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
  int *temp_idx = (int*)malloc(n * sizeof(int));
  for(size_t k = 0; k < n; k++) temp_idx[k] = (int)k;
  srand(0);
  for(int i_rand = 0; i_rand < nrand; i_rand++) {
    shuffle(temp_idx, (int)n);
    memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
  }
  free(temp_idx);

  size_t total_elements = p * m;
  int nths = 1;
  #ifdef _OPENMP
    #pragma omp parallel
    { 
       // FIXED: Multiline pragma to satisfy Clang expression parser
       #pragma omp single
       nths = omp_get_num_threads(); 
    }
  #endif

  double *tsum = (double*)calloc(nths * total_elements, sizeof(double));
  double *tsum_sq = (double*)calloc(nths * total_elements, sizeof(double));
  double *tcount = (double*)calloc(nths * total_elements, sizeof(double));

  #pragma omp parallel
  {
    int tid = 0;
    #ifdef _OPENMP
      tid = omp_get_thread_num();
    #endif
    double *my_s = &tsum[tid * total_elements], *my_ss = &tsum_sq[tid * total_elements], *my_c = &tcount[tid * total_elements];
    gsl_matrix *Tt_p = gsl_matrix_alloc(n, p);
    gsl_matrix *br = gsl_matrix_alloc(p, m);

    #pragma omp for schedule(dynamic, 16)
    for(int ir = 0; ir < nrand; ir++) {
      int *p_idx = &perm_table[(size_t)ir * n];
      for(size_t i = 0; i < n; i++) memcpy(Tt_p->data + ((size_t)p_idx[i] * p), Tt_orig->data + (i * p), p * sizeof(double));
      gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Tt_p, Yt, 0.0, br);
      for(size_t k = 0; k < total_elements; k++) {
        double val = br->data[k];
        my_s[k] += val; my_ss[k] += val * val;
        if(fabs(val) >= fabs(beta_vec[k])) my_c[k] += 1.0;
      }
    }
    gsl_matrix_free(Tt_p); gsl_matrix_free(br);
  }

  #pragma omp parallel for schedule(static)
  for(size_t i = 0; i < total_elements; i++) {
    double s = 0, ss = 0, c = 0;
    for(int t = 0; t < nths; t++) { s += tsum[t*total_elements+i]; ss += tsum_sq[t*total_elements+i]; c += tcount[t*total_elements+i]; }
    zscore_vec[i] = s; se_vec[i] = ss; pvalue_vec[i] = c;
  }

  double inv_nr = 1.0 / nrand;
  #pragma omp parallel for schedule(static)
  for(size_t i = 0; i < total_elements; i++) {
    pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
    double mean = zscore_vec[i] * inv_nr;
    double sd = sqrt(fmax(0, (se_vec[i] * inv_nr) - (mean * mean)));
    se_vec[i] = sd;
    zscore_vec[i] = (sd > 1e-12) ? (beta_vec[i] - mean) / sd : 0.0;
  }

  free(tsum); free(tsum_sq); free(tcount); free(perm_table); gsl_matrix_free(Tt_orig);
  gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_partial_free(Xt); 
  gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta);
}

/* =========================================================
   INTERFACE WRAPPERS
   ========================================================= */

SEXP make_res_list(SEXP b, SEXP s, SEXP z, SEXP p) {
  SEXP res = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, b); SET_VECTOR_ELT(res, 1, s); SET_VECTOR_ELT(res, 2, z); SET_VECTOR_ELT(res, 3, p);
  SEXP names = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(names, 0, mkChar("beta")); SET_STRING_ELT(names, 1, mkChar("se"));
  SET_STRING_ELT(names, 2, mkChar("zscore")); SET_STRING_ELT(names, 3, mkChar("pvalue"));
  setAttrib(res, R_NamesSymbol, names);
  UNPROTECT(2); return res;
}

SEXP ridgeReg_old_interface(SEXP X, SEXP Y, SEXP L, SEXP R) {
  size_t n = INTEGER(getAttrib(X, R_DimSymbol))[0], p = INTEGER(getAttrib(X, R_DimSymbol))[1], m = INTEGER(getAttrib(Y, R_DimSymbol))[1];
  SEXP b = PROTECT(allocVector(REALSXP, p*m)), s = PROTECT(allocVector(REALSXP, p*m)), z = PROTECT(allocVector(REALSXP, p*m)), pv = PROTECT(allocVector(REALSXP, p*m));
  ridgeReg_old_core(REAL(X), REAL(Y), n, p, m, asReal(L), asInteger(R), REAL(b), REAL(s), REAL(z), REAL(pv));
  SEXP res = make_res_list(b, s, z, pv); UNPROTECT(4); return res;
}

SEXP ridgeRegTperm_old_interface(SEXP X, SEXP Y, SEXP L, SEXP R) {
  size_t n = INTEGER(getAttrib(X, R_DimSymbol))[0], p = INTEGER(getAttrib(X, R_DimSymbol))[1], m = INTEGER(getAttrib(Y, R_DimSymbol))[1];
  SEXP b = PROTECT(allocVector(REALSXP, p*m)), s = PROTECT(allocVector(REALSXP, p*m)), z = PROTECT(allocVector(REALSXP, p*m)), pv = PROTECT(allocVector(REALSXP, p*m));
  ridgeRegTperm_old_core(REAL(X), REAL(Y), n, p, m, asReal(L), asInteger(R), REAL(b), REAL(s), REAL(z), REAL(pv));
  SEXP res = make_res_list(b, s, z, pv); UNPROTECT(4); return res;
}

SEXP ridgeRegFast_interface(SEXP X, SEXP Y, SEXP L, SEXP R, SEXP C) {
  size_t n = INTEGER(getAttrib(X, R_DimSymbol))[0], p = INTEGER(getAttrib(X, R_DimSymbol))[1], m = INTEGER(getAttrib(Y, R_DimSymbol))[1];
  SEXP b = PROTECT(allocVector(REALSXP, p*m)), s = PROTECT(allocVector(REALSXP, p*m)), z = PROTECT(allocVector(REALSXP, p*m)), pv = PROTECT(allocVector(REALSXP, p*m));
  ridgeRegFast_core(REAL(X), REAL(Y), n, p, m, asReal(L), asInteger(R), REAL(b), REAL(s), REAL(z), REAL(pv), asInteger(C));
  SEXP res = make_res_list(b, s, z, pv); UNPROTECT(4); return res;
}

SEXP ridgeRegTperm_interface(SEXP X, SEXP Y, SEXP L, SEXP R, SEXP C) {
  size_t n = INTEGER(getAttrib(X, R_DimSymbol))[0], p = INTEGER(getAttrib(X, R_DimSymbol))[1], m = INTEGER(getAttrib(Y, R_DimSymbol))[1];
  SEXP b = PROTECT(allocVector(REALSXP, p*m)), s = PROTECT(allocVector(REALSXP, p*m)), z = PROTECT(allocVector(REALSXP, p*m)), pv = PROTECT(allocVector(REALSXP, p*m));
  ridgeRegTperm_core(REAL(X), REAL(Y), n, p, m, asReal(L), asInteger(R), REAL(b), REAL(s), REAL(z), REAL(pv), asInteger(C));
  SEXP res = make_res_list(b, s, z, pv); UNPROTECT(4); return res;
}

static const R_CallMethodDef callMethods[] = {
  {"ridgeReg_old_interface",      (DL_FUNC) &ridgeReg_old_interface,      4},
  {"ridgeRegTperm_old_interface", (DL_FUNC) &ridgeRegTperm_old_interface, 4},
  {"ridgeRegFast_interface",      (DL_FUNC) &ridgeRegFast_interface,      5},
  {"ridgeRegTperm_interface",     (DL_FUNC) &ridgeRegTperm_interface,     5},
  {NULL, NULL, 0}
};

void R_init_RidgeR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
