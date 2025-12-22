#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

/* =========================================================
   HELPER FUNCTIONS
   ========================================================= */

// Custom GSL matrix wrapper for R-allocated memory
gsl_matrix *RVectorObject_to_gsl_matrix(double *vec, size_t nr, size_t nc) {
    gsl_block *b = (gsl_block*)malloc(sizeof(gsl_block));
    gsl_matrix *r = (gsl_matrix*)malloc(sizeof(gsl_matrix));
    r->size1 = nr; r->tda = r->size2 = nc; r->owner = 1; 
    b->data = r->data = vec; r->block = b; b->size = nr * nc;
    return r;
}

void gsl_matrix_partial_free(gsl_matrix *x) {
    if(x) { if(x->block) free(x->block); free(x); }
}

void shuffle(int array[], const int n) {
    int i, j, t;
    for (i = 0; i < n - 1; i++) {
        j = i + rand() / (RAND_MAX / (n - i) + 1);
        t = array[j]; array[j] = array[i]; array[i] = t;
    }
}

/* =========================================================
   VERSION 1: ORIGINAL (Single-threaded, Y row permutation)
   ========================================================= */
void ridgeReg(double *X_vec, double *Y_vec, int *n_pt, int *p_pt, int *m_pt,
              double *lambda_pt, double *nrand_pt,
              double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec) {
    int n = *n_pt, p = *p_pt, m = *m_pt, nrand = (int)*nrand_pt;
    double lambda = *lambda_pt;

    gsl_matrix *X = RVectorObject_to_gsl_matrix(X_vec, n, p);
    gsl_matrix *Y = RVectorObject_to_gsl_matrix(Y_vec, n, m);
    gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);
    
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y_rand = gsl_matrix_alloc(n, m), *beta_rand = gsl_matrix_alloc(p, m);
    gsl_matrix *aver = gsl_matrix_alloc(p, m), *aver_sq = gsl_matrix_alloc(p, m);

    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);
    gsl_linalg_cholesky_decomp(I); gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

    int *array_index = (int*)malloc(n * sizeof(int));
    for(int i=0; i<n; i++) array_index[i] = i;
    srand(0);

    gsl_matrix_set_zero(aver); gsl_matrix_set_zero(aver_sq);
    double *pv = pvalue_vec; 
    size_t total = (size_t)p * m;
    for(size_t i=0; i<total; i++) pv[i] = 0.0;

    for(int i=0; i<nrand; i++) {
        shuffle(array_index, n);
        for(int j=0; j<n; j++){
            gsl_vector_const_view t_row = gsl_matrix_const_row(Y, array_index[j]);
            gsl_matrix_set_row(Y_rand, j, &t_row.vector);
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);
        for(size_t j=0; j<total; j++) {
            if(fabs(beta_rand->data[j]) >= fabs(beta->data[j])) pv[j]++;
        }
        gsl_matrix_add(aver, beta_rand);
        gsl_matrix_mul_elements(beta_rand, beta_rand);
        gsl_matrix_add(aver_sq, beta_rand);
    }

    double inv_n = 1.0 / nrand;
    for(size_t i=0; i<total; i++) {
        pvalue_vec[i] = (pv[i] + 1.0) / (nrand + 1.0);
        double m_r = aver->data[i] * inv_n;
        double v_r = (aver_sq->data[i] * inv_n) - (m_r * m_r);
        double s_r = sqrt(fmax(0.0, v_r));
        se_vec[i] = s_r;
        zscore_vec[i] = (s_r > 1e-12) ? (beta_vec[i] - m_r) / s_r : 0.0;
    }

    free(array_index); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_rand);
    gsl_matrix_free(beta_rand); gsl_matrix_free(aver); gsl_matrix_free(aver_sq);
    gsl_matrix_partial_free(X); gsl_matrix_partial_free(Y); gsl_matrix_partial_free(beta);
}

/* =========================================================
   VERSION 3: T COLUMN PERMUTATION (Multi-threaded Core)
   ========================================================= */
void ridgeRegTperm_core(double *X_ptr, double *Y_ptr, size_t n, size_t p, size_t m,
                        double lambda, int nrand, double *beta_vec, double *se_vec, 
                        double *zscore_vec, double *pvalue_vec, int num_threads) {
    #ifdef _OPENMP
    if (num_threads > 0) omp_set_num_threads(num_threads);
    #endif

    gsl_matrix *Xt = RVectorObject_to_gsl_matrix(X_ptr, p, n);
    gsl_matrix *Yt = RVectorObject_to_gsl_matrix(Y_ptr, m, n);
    gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);
    gsl_matrix *I_mat = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I_mat);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
    gsl_linalg_cholesky_decomp(I_mat); gsl_linalg_cholesky_invert(I_mat);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);
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

    R_xlen_t total_elements = (R_xlen_t)p * m;
    int actual_threads = 1;
    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        actual_threads = omp_get_num_threads();
    }
    #endif

    double *th_sum = (double*)calloc((size_t)actual_threads * total_elements, sizeof(double));
    double *th_sq  = (double*)calloc((size_t)actual_threads * total_elements, sizeof(double));
    double *th_cnt = (double*)calloc((size_t)actual_threads * total_elements, sizeof(double));

    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        double *m_s = &th_sum[(size_t)tid * total_elements];
        double *m_q = &th_sq[(size_t)tid * total_elements];
        double *m_c = &th_cnt[(size_t)tid * total_elements];

        gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p), *br = gsl_matrix_alloc(p, m);
        size_t rsz = p * sizeof(double);

        #pragma omp for schedule(dynamic, 16)
        for(int i_rand = 0; i_rand < nrand; i_rand++) {
            int *p_idx = &perm_table[(size_t)i_rand * n];
            for(size_t i = 0; i < n; i++) 
                memcpy(Tt_perm->data + ((size_t)p_idx[i] * p), Tt_orig->data + (i * p), rsz);
            
            gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Tt_perm, Yt, 0.0, br);

            for(size_t k = 0; k < total_elements; k++) {
                double b_r = br->data[k];
                m_s[k] += b_r; m_q[k] += b_r * b_r;
                if(fabs(b_r) >= fabs(beta_vec[k])) m_c[k] += 1.0;
            }
        }
        gsl_matrix_free(Tt_perm); gsl_matrix_free(br);
    }

    double inv_n = 1.0 / nrand;
    #pragma omp parallel for schedule(static)
    for(R_xlen_t i = 0; i < total_elements; i++) {
        double s=0, q=0, c=0;
        for(int t=0; t<actual_threads; t++) {
            size_t off = (size_t)t * total_elements + i;
            s += th_sum[off]; q += th_sq[off]; c += th_cnt[off];
        }
        pvalue_vec[i] = (c + 1.0) / (nrand + 1.0);
        double mean = s * inv_n;
        double var  = (q * inv_n) - (mean * mean);
        double std  = sqrt(fmax(0.0, var));
        se_vec[i] = std;
        zscore_vec[i] = (std > 1e-12) ? (beta_vec[i] - mean) / std : 0.0;
    }

    free(th_sum); free(th_sq); free(th_cnt); free(perm_table);
    gsl_matrix_free(Tt_orig); gsl_matrix_free(I_mat); gsl_matrix_free(T);
    gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta);
}

SEXP ridgeRegTperm_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp, SEXP ncores_sexp) {
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    size_t m = (size_t)INTEGER(getAttrib(Y_sexp, R_DimSymbol))[1];
    R_xlen_t total_len = (R_xlen_t)p * m;

    SEXP b_s = PROTECT(allocVector(REALSXP, total_len)), s_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP z_s = PROTECT(allocVector(REALSXP, total_len)), p_s = PROTECT(allocVector(REALSXP, total_len));

    ridgeRegTperm_core(REAL(X_sexp), REAL(Y_sexp), n, p, m, asReal(lambda_sexp), 
                       asInteger(nrand_sexp), REAL(b_s), REAL(s_s), REAL(z_s), REAL(p_s), 
                       asInteger(ncores_sexp));

    SEXP res = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(res, 0, b_s); SET_VECTOR_ELT(res, 1, s_s);
    SET_VECTOR_ELT(res, 2, z_s); SET_VECTOR_ELT(res, 3, p_s);
    SEXP nms = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, mkChar("beta")); SET_STRING_ELT(nms, 1, mkChar("se"));
    SET_STRING_ELT(nms, 2, mkChar("zscore")); SET_STRING_ELT(nms, 3, mkChar("pvalue"));
    setAttrib(res, R_NamesSymbol, nms);

    UNPROTECT(6); return res;
}

/* =========================================================
   REGISTRATION
   ========================================================= */
static const R_CMethodDef cMethods[] = {
    {"ridgeReg", (DL_FUNC) &ridgeReg, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef callMethods[] = {
    {"ridgeRegTperm_interface", (DL_FUNC) &ridgeRegTperm_interface, 5},
    {NULL, NULL, 0}
};

void R_init_RidgeR(DllInfo *dll) {
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
