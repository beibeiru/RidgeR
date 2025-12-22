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

void shuffle(int array[], const int n) {
    int i, j, t;
    for (i = 0; i < n - 1; i++) {
        j = i + rand() / (RAND_MAX / (n - i) + 1);
        t = array[j]; array[j] = array[i]; array[i] = t;
    }
}

int* generate_perm_table(int n, int nrand) {
    int *table = (int*)malloc((size_t)nrand * (size_t)n * sizeof(int));
    int *base_idx = (int*)malloc((size_t)n * sizeof(int));
    for(int i=0; i<n; i++) base_idx[i] = i;
    srand(0); 
    for(int i=0; i<nrand; i++) {
        shuffle(base_idx, n);
        memcpy(table + ((size_t)i * n), base_idx, n * sizeof(int));
    }
    free(base_idx);
    return table;
}

gsl_matrix *R_to_gsl_matrix(double *vec, size_t nr, size_t nc) {
    gsl_block *b = (gsl_block*)malloc(sizeof(gsl_block));
    gsl_matrix *r = (gsl_matrix*)malloc(sizeof(gsl_matrix));
    r->size1 = nr; r->tda = r->size2 = nc; r->owner = 1; 
    b->data = r->data = vec; r->block = b; b->size = nr * nc;
    return r;
}

void gsl_matrix_partial_free(gsl_matrix *x) {
    if(x) { if(x->block) free(x->block); free(x); }
}

SEXP create_res_list(SEXP beta, SEXP se, SEXP zscore, SEXP pvalue) {
    SEXP res_list = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(res_list, 0, beta); SET_VECTOR_ELT(res_list, 1, se);
    SET_VECTOR_ELT(res_list, 2, zscore); SET_VECTOR_ELT(res_list, 3, pvalue);
    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("beta")); SET_STRING_ELT(names, 1, mkChar("se"));
    SET_STRING_ELT(names, 2, mkChar("zscore")); SET_STRING_ELT(names, 3, mkChar("pvalue"));
    setAttrib(res_list, R_NamesSymbol, names);
    UNPROTECT(2); return res_list;
}

/* =========================================================
   VERSION 0: LEGACY .C (ridgeReg)
   ========================================================= */
void ridgeReg(double *X_vec, double *Y_vec, int *n_pt, int *p_pt, int *m_pt,
              double *lambda_pt, double *nrand_pt,
              double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec) {
    size_t n = (size_t)*n_pt, p = (size_t)*p_pt, m = (size_t)*m_pt;
    int nrand = (int)*nrand_pt; double lambda = *lambda_pt;

    gsl_matrix *Xt = R_to_gsl_matrix(X_vec, p, n), *Yt = R_to_gsl_matrix(Y_vec, m, n);
    gsl_matrix *beta = R_to_gsl_matrix(beta_vec, p, m);
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Yr = gsl_matrix_alloc(n, m), *Yp = gsl_matrix_alloc(n, m), *br = gsl_matrix_alloc(p, m);
    gsl_matrix *aver = gsl_matrix_alloc(p, m), *aver_sq = gsl_matrix_alloc(p, m);

    gsl_matrix_set_identity(I); gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I);
    gsl_linalg_cholesky_decomp(I); gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, Xt, 0.0, T);
    gsl_matrix_transpose_memcpy(Yr, Yt);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, beta);

    int *p_tab = generate_perm_table((int)n, nrand);
    gsl_matrix_set_zero(aver); gsl_matrix_set_zero(aver_sq);
    size_t total = p * m; for(size_t i=0; i<total; i++) pvalue_vec[i] = 0;

    for(int i=0; i<nrand; i++) {
        int *cur = p_tab + (i * n);
        for(size_t g=0; g<n; g++) memcpy(Yp->data + (g * m), Yr->data + ((size_t)cur[g] * m), m * sizeof(double));
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yp, 0.0, br);
        for(size_t j=0; j<total; j++) if(fabs(br->data[j]) >= fabs(beta->data[j])) pvalue_vec[j]++;
        gsl_matrix_add(aver, br); gsl_matrix_mul_elements(br, br); gsl_matrix_add(aver_sq, br);
    }

    double inv_n = 1.0 / nrand;
    for(size_t i=0; i<total; i++) {
        pvalue_vec[i] = (pvalue_vec[i] + 1.0) / (nrand + 1.0);
        double mean = aver->data[i] * inv_n; double var = (aver_sq->data[i] * inv_n) - (mean * mean);
        double std = sqrt(fmax(0.0, var));
        se_vec[i] = std; zscore_vec[i] = (std > 1e-12) ? (beta_vec[i] - mean) / std : 0.0;
    }
    free(p_tab); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Yr); gsl_matrix_free(Yp);
    gsl_matrix_free(br); gsl_matrix_free(aver); gsl_matrix_free(aver_sq);
    gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta);
}

/* =========================================================
   VERSION 1: OLD (.Call Single-threaded Y-Perm)
   ========================================================= */
SEXP ridgeReg_old_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s) {
    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    size_t m = (size_t)INTEGER(getAttrib(Y_s, R_DimSymbol))[1];
    int nrand = asInteger(nrand_s); double lambda = asReal(lambda_s);
    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len)), s_s = PROTECT(allocVector(REALSXP, len));
    SEXP z_s = PROTECT(allocVector(REALSXP, len)), p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);
    memset(pv, 0, len * sizeof(double)); memset(sv, 0, len * sizeof(double)); memset(zv, 0, len * sizeof(double));

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_s), p, n), Yt = gsl_matrix_view_array(REAL(Y_s), m, n);
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n), *Yr = gsl_matrix_alloc(n, m);
    gsl_matrix *Yp = gsl_matrix_alloc(n, m), *br = gsl_matrix_alloc(p, m);

    gsl_matrix_set_identity(I); gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I); gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);
    gsl_matrix_transpose_memcpy(Yr, &Yt.matrix);
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, &bo.matrix);

    int *p_tab = generate_perm_table((int)n, nrand);
    for(int i=0; i<nrand; i++) {
        int *cur = p_tab + (i * n);
        for(size_t g=0; g<n; g++) memcpy(Yp->data + (g * m), Yr->data + ((size_t)cur[g] * m), m * sizeof(double));
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yp, 0.0, br);
        for(R_xlen_t k=0; k<len; k++) {
            if(fabs(br->data[k]) >= fabs(bv[k])) pv[k]++;
            zv[k] += br->data[k]; sv[k] += br->data[k] * br->data[k];
        }
    }
    double inv_n = 1.0 / nrand;
    for(R_xlen_t i=0; i<len; i++) {
        double mean = zv[i] * inv_n; double var = (sv[i] * inv_n) - (mean * mean);
        sv[i] = sqrt(fmax(0.0, var)); zv[i] = (sv[i] > 1e-12) ? (bv[i] - mean) / sv[i] : 0.0;
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
    }
    free(p_tab); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Yr); gsl_matrix_free(Yp); gsl_matrix_free(br);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s); UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 2: OLD1 (.Call Single-threaded T-Perm Scatter)
   ========================================================= */
SEXP ridgeRegTperm_old_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s) {
    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    size_t m = (size_t)INTEGER(getAttrib(Y_s, R_DimSymbol))[1];
    int nrand = asInteger(nrand_s); double lambda = asReal(lambda_s);
    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len)), s_s = PROTECT(allocVector(REALSXP, len));
    SEXP z_s = PROTECT(allocVector(REALSXP, len)), p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);
    memset(pv, 0, len * sizeof(double)); memset(sv, 0, len * sizeof(double)); memset(zv, 0, len * sizeof(double));

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_s), p, n), Yt = gsl_matrix_view_array(REAL(Y_s), m, n);
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n), *Yr = gsl_matrix_alloc(n, m);
    gsl_matrix_set_identity(I); gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I); gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);
    gsl_matrix_transpose_memcpy(Yr, &Yt.matrix);
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, &bo.matrix);

    gsl_matrix *Tt_o = gsl_matrix_alloc(n, p), *Ttp = gsl_matrix_alloc(n, p), *br = gsl_matrix_alloc(p, m);
    gsl_matrix_transpose_memcpy(Tt_o, T);
    int *p_tab = generate_perm_table((int)n, nrand);

    for(int i_r=0; i_r<nrand; i_r++) {
        int *pidx = p_tab + (i_r * n);
        for(size_t i=0; i<n; i++) memcpy(Ttp->data + ((size_t)pidx[i] * p), Tt_o->data + (i * p), p * sizeof(double));
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Ttp, Yr, 0.0, br);
        for(R_xlen_t k=0; k<len; k++) {
            if(fabs(br->data[k]) >= fabs(bv[k])) pv[k]++;
            zv[k] += br->data[k]; sv[k] += br->data[k] * br->data[k];
        }
    }
    double inv_n = 1.0 / nrand;
    for(R_xlen_t i=0; i<len; i++) {
        double mean = zv[i] * inv_n; double var = (sv[i] * inv_n) - (mean * mean);
        sv[i] = sqrt(fmax(0.0, var)); zv[i] = (sv[i] > 1e-12) ? (bv[i] - mean) / sv[i] : 0.0;
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
    }
    free(p_tab); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Yr); gsl_matrix_free(Tt_o); gsl_matrix_free(Ttp); gsl_matrix_free(br);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s); UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 3: NEW (.Call Multi-threaded Y-Perm)
   ========================================================= */
SEXP ridgeRegFast_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s, SEXP ncores_s) {
    int ncores = asInteger(ncores_s);
    #ifdef _OPENMP
      if (ncores > 0) omp_set_num_threads(ncores);
    #endif
    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1], m = (size_t)INTEGER(getAttrib(Y_s, R_DimSymbol))[1];
    int nrand = asInteger(nrand_s); double lambda = asReal(lambda_s);
    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len)), s_s = PROTECT(allocVector(REALSXP, len)), z_s = PROTECT(allocVector(REALSXP, len)), p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);
    memset(pv, 0, len * sizeof(double)); memset(sv, 0, len * sizeof(double)); memset(zv, 0, len * sizeof(double));

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_s), p, n), Yt = gsl_matrix_view_array(REAL(Y_s), m, n);
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n), *Yr = gsl_matrix_alloc(n, m);
    gsl_matrix_set_identity(I); gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I); gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);
    gsl_matrix_transpose_memcpy(Yr, &Yt.matrix);
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, &bo.matrix);

    int *p_tab = generate_perm_table((int)n, nrand);
    int actual_threads = 1;

    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        {
            actual_threads = omp_get_num_threads();
        }
    }
    #endif

    double *th_s = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *th_q = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *th_c = (double*)calloc((size_t)actual_threads * len, sizeof(double));

    #pragma omp parallel
    {
        int tid = 0; 
        #ifdef _OPENMP
        tid = omp_get_thread_num(); 
        #endif

        double *ms = &th_s[(size_t)tid * len];
        double *mq = &th_q[(size_t)tid * len];
        double *mc = &th_c[(size_t)tid * len];

        gsl_matrix *Yp = gsl_matrix_alloc(n, m), *br = gsl_matrix_alloc(p, m);
        
        #pragma omp for schedule(dynamic)
        for(int i = 0; i < nrand; i++) {
            int *pidx = p_tab + ((size_t)i * n);
            for(size_t g = 0; g < n; g++) memcpy(Yp->data + (g * m), Yr->data + ((size_t)pidx[g] * m), m * sizeof(double));
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yp, 0.0, br);
            for(R_xlen_t k=0; k<len; k++) {
                if(fabs(br->data[k]) >= fabs(bv[k])) mc[k]++;
                ms[k] += br->data[k]; mq[k] += br->data[k] * br->data[k];
            }
        }
        gsl_matrix_free(Yp); gsl_matrix_free(br);
    }

    for(int t=0; t<actual_threads; t++) {
        for(R_xlen_t i=0; i<len; i++) { 
            zv[i] += th_s[(size_t)t*len+i]; 
            sv[i] += th_q[(size_t)t*len+i]; 
            pv[i] += th_c[(size_t)t*len+i]; 
        }
    }

    double inv_n = 1.0 / nrand;
    for(R_xlen_t i=0; i<len; i++) {
        double mean = zv[i] * inv_n; double var = (sv[i] * inv_n) - (mean * mean);
        sv[i] = sqrt(fmax(0.0, var)); zv[i] = (sv[i] > 1e-12) ? (bv[i] - mean) / sv[i] : 0.0;
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
    }

    free(th_s); free(th_q); free(th_c); free(p_tab); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Yr);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s); UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 4: NEW2 (.Call Multi-threaded T-Perm Scatter)
   ========================================================= */
SEXP ridgeRegTperm_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s, SEXP ncores_s) {
    int ncores = asInteger(ncores_s);
    #ifdef _OPENMP
      if (ncores > 0) omp_set_num_threads(ncores);
    #endif
    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1], m = (size_t)INTEGER(getAttrib(Y_s, R_DimSymbol))[1];
    int nrand = asInteger(nrand_s); double lambda = asReal(lambda_s);
    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len)), s_s = PROTECT(allocVector(REALSXP, len)), z_s = PROTECT(allocVector(REALSXP, len)), p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);
    memset(pv, 0, len * sizeof(double)); memset(sv, 0, len * sizeof(double)); memset(zv, 0, len * sizeof(double));

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_s), p, n), Yt = gsl_matrix_view_array(REAL(Y_s), m, n);
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n), *Yr = gsl_matrix_alloc(n, m);
    gsl_matrix_set_identity(I); gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I); gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);
    gsl_matrix_transpose_memcpy(Yr, &Yt.matrix);
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Yr, 0.0, &bo.matrix);

    gsl_matrix *Tt_o = gsl_matrix_alloc(n, p); gsl_matrix_transpose_memcpy(Tt_o, T);
    int *p_tab = generate_perm_table((int)n, nrand);
    int actual_threads = 1;

    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        {
            actual_threads = omp_get_num_threads();
        }
    }
    #endif

    double *th_s = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *th_q = (double*)calloc((size_t)actual_threads * len, sizeof(double));
    double *th_c = (double*)calloc((size_t)actual_threads * len, sizeof(double));

    #pragma omp parallel
    {
        int tid = 0; 
        #ifdef _OPENMP
        tid = omp_get_thread_num(); 
        #endif

        double *ms = &th_s[(size_t)tid * len];
        double *mq = &th_q[(size_t)tid * len];
        double *mc = &th_c[(size_t)tid * len];

        gsl_matrix *Ttp = gsl_matrix_alloc(n, p), *br = gsl_matrix_alloc(p, m);
        
        #pragma omp for schedule(dynamic, 16)
        for(int i_r = 0; i_r < nrand; i_r++) {
            int *pidx = p_tab + ((size_t)i_r * n);
            for(size_t i = 0; i < n; i++) memcpy(Ttp->data + ((size_t)pidx[i] * p), Tt_o->data + (i * p), p * sizeof(double));
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Ttp, Yr, 0.0, br);
            for(R_xlen_t k=0; k<len; k++) {
                if(fabs(br->data[k]) >= fabs(bv[k])) mc[k]++;
                ms[k] += br->data[k]; mq[k] += br->data[k] * br->data[k];
            }
        }
        gsl_matrix_free(Ttp); gsl_matrix_free(br);
    }

    for(int t=0; t<actual_threads; t++) {
        for(R_xlen_t i=0; i<len; i++) { 
            zv[i] += th_s[(size_t)t*len+i]; 
            sv[i] += th_q[(size_t)t*len+i]; 
            pv[i] += th_c[(size_t)t*len+i]; 
        }
    }

    double inv_n = 1.0 / nrand;
    for(R_xlen_t i = 0; i < len; i++) {
        double mean = zv[i] * inv_n; double var = (sv[i] * inv_n) - (mean * mean);
        sv[i] = sqrt(fmax(0.0, var)); zv[i] = (sv[i] > 1e-12) ? (bv[i] - mean) / sv[i] : 0.0;
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
    }
    free(th_s); free(th_q); free(th_c); free(p_tab); gsl_matrix_free(Tt_o); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Yr);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s); UNPROTECT(4); return res;
}

static const R_CMethodDef cMethods[] = { {"ridgeReg", (DL_FUNC) &ridgeReg, 11}, {NULL, NULL, 0} };
static const R_CallMethodDef callMethods[] = {
    {"ridgeReg_old_interface", (DL_FUNC) &ridgeReg_old_interface, 4},
    {"ridgeRegTperm_old_interface", (DL_FUNC) &ridgeRegTperm_old_interface, 4},
    {"ridgeRegFast_interface", (DL_FUNC) &ridgeRegFast_interface, 5},
    {"ridgeRegTperm_interface", (DL_FUNC) &ridgeRegTperm_interface, 5},
    {NULL, NULL, 0}
};
void R_init_RidgeR(DllInfo *dll) { R_registerRoutines(dll, cMethods, callMethods, NULL, NULL); R_useDynamicSymbols(dll, FALSE); }
