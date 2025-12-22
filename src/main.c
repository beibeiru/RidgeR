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

// Tuning parameters for the Fast version
#define SAMPLE_STRIP_SIZE 64
#define PERM_BATCH_SIZE 64

/* =========================================================
   HELPER FUNCTIONS
   ========================================================= */

void shuffle(int array[], const int n) {
    int i, j, t;
    for (i = 0; i < n - 1; i++) {
        j = i + rand() / (RAND_MAX / (n - i) + 1);
        t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

// Wraps the four result vectors into a named R list
SEXP create_res_list(SEXP beta, SEXP se, SEXP zscore, SEXP pvalue) {
    SEXP res_list = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(res_list, 0, beta);
    SET_VECTOR_ELT(res_list, 1, se);
    SET_VECTOR_ELT(res_list, 2, zscore);
    SET_VECTOR_ELT(res_list, 3, pvalue);
    
    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("beta"));
    SET_STRING_ELT(names, 1, mkChar("se"));
    SET_STRING_ELT(names, 2, mkChar("zscore"));
    SET_STRING_ELT(names, 3, mkChar("pvalue"));
    setAttrib(res_list, R_NamesSymbol, names);
    
    UNPROTECT(2);
    return res_list;
}

/* =========================================================
   VERSION 1: OLD (Single-threaded, Y Row Permutation)
   ========================================================= */
SEXP ridgeReg_old_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s) {
    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_s, R_DimSymbol);
    size_t m = (size_t)INTEGER(y_dim)[1];
    int nrand = asInteger(nrand_s);
    double lambda = asReal(lambda_s);

    R_xlen_t total_len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP s_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP z_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP p_s = PROTECT(allocVector(REALSXP, total_len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_s), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_s), m, n);
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    
    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);
    
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, &Yt.matrix, 0.0, &bo.matrix);

    int *idx = (int*)malloc(n * sizeof(int));
    for(size_t i=0; i<n; i++) idx[i] = (int)i;
    srand(0);
    for(R_xlen_t i=0; i<total_len; i++) { pv[i]=0; sv[i]=0; zv[i]=0; }

    gsl_matrix *Y_perm = gsl_matrix_alloc(n, m);
    gsl_matrix *br = gsl_matrix_alloc(p, m);

    for(int i=0; i<nrand; i++) {
        shuffle(idx, (int)n);
        for(size_t g=0; g<n; g++) {
            for(size_t s=0; s<m; s++) {
                gsl_matrix_set(Y_perm, g, s, gsl_matrix_get(&Yt.matrix, s, idx[g]));
            }
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_perm, 0.0, br);
        for(R_xlen_t k=0; k < total_len; k++) {
            if(fabs(br->data[k]) >= fabs(bv[k])) pv[k] += 1.0;
            zv[k] += br->data[k]; sv[k] += br->data[k] * br->data[k];
        }
    }

    double inv_n = 1.0 / nrand;
    for(R_xlen_t i=0; i < total_len; i++) {
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
        double mr = zv[i] * inv_n;
        double vr = (sv[i] * inv_n) - (mr * mr);
        double sr = sqrt(fmax(0, vr));
        sv[i] = sr;
        zv[i] = (sr > 1e-12) ? (bv[i] - mr) / sr : 0.0;
    }

    free(idx); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_perm); gsl_matrix_free(br);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 1.5: OLD2 (Single-threaded, T Column Permutation)
   ========================================================= */
SEXP ridgeRegTperm_old_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s) {
    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_s, R_DimSymbol);
    size_t m = (size_t)INTEGER(y_dim)[1];
    int nrand = asInteger(nrand_s);
    double lambda = asReal(lambda_s);

    R_xlen_t total_len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP s_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP z_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP p_s = PROTECT(allocVector(REALSXP, total_len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_s), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_s), m, n);
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    
    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);
    
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, &Yt.matrix, 0.0, &bo.matrix);

    gsl_matrix *Tt_o = gsl_matrix_alloc(n, p);
    gsl_matrix_transpose_memcpy(Tt_o, T);
    gsl_matrix *Tt_p = gsl_matrix_alloc(n, p);
    gsl_matrix *br = gsl_matrix_alloc(p, m);

    int *p_tab = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *tmp = (int*)malloc(n * sizeof(int));
    for(size_t k=0; k<n; k++) tmp[k] = (int)k;
    srand(0);
    for(int i=0; i<nrand; i++) { shuffle(tmp, (int)n); memcpy(&p_tab[(size_t)i * n], tmp, n * sizeof(int)); }
    free(tmp);

    for(R_xlen_t i=0; i<total_len; i++) { pv[i]=0; sv[i]=0; zv[i]=0; }
    size_t row_sz = p * sizeof(double);

    for(int i_r=0; i_r<nrand; i_r++) {
        int *p_idx = &p_tab[(size_t)i_r * n];
        for(size_t i=0; i<n; i++) memcpy(Tt_p->data + (p_idx[i] * p), Tt_o->data + (i * p), row_sz);
        gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Tt_p, &Yt.matrix, 0.0, br);
        for(R_xlen_t k=0; k<total_len; k++) {
            if(fabs(br->data[k]) >= fabs(bv[k])) pv[k] += 1.0;
            zv[k] += br->data[k]; sv[k] += br->data[k] * br->data[k];
        }
    }

    double inv_n = 1.0 / nrand;
    for(R_xlen_t i=0; i<total_len; i++) {
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
        double mr = zv[i] * inv_n;
        double vr = (sv[i] * inv_n) - (mr * mr);
        double sr = sqrt(fmax(0, vr));
        sv[i] = sr; zv[i] = (sr > 1e-12) ? (bv[i] - mr) / sr : 0.0;
    }

    free(p_tab); gsl_matrix_free(Tt_o); gsl_matrix_free(Tt_p); gsl_matrix_free(br);
    gsl_matrix_free(I); gsl_matrix_free(T);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 2: NEW (Multi-threaded, Y Row Permutation)
   ========================================================= */
SEXP ridgeRegFast_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s, SEXP ncores_s) {
    int ncores = asInteger(ncores_s);
    #ifdef _OPENMP
    if (ncores > 0) omp_set_num_threads(ncores);
    #endif

    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_s, R_DimSymbol);
    size_t m = (size_t)INTEGER(y_dim)[1];
    double lambda = asReal(lambda_s);
    int nrand = asInteger(nrand_s);

    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len));
    SEXP s_s = PROTECT(allocVector(REALSXP, len));
    SEXP z_s = PROTECT(allocVector(REALSXP, len));
    SEXP p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_s), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_s), m, n);
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);
    
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, &Yt.matrix, 0.0, &bo.matrix);

    int *p_tab = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *tmp = (int*)malloc(n * sizeof(int));
    for(size_t k=0; k<n; k++) tmp[k] = (int)k;
    srand(0);
    for(int i=0; i<nrand; i++) { shuffle(tmp, (int)n); memcpy(&p_tab[(size_t)i * n], tmp, n * sizeof(int)); }
    free(tmp);

    #pragma omp parallel for schedule(static)
    for(R_xlen_t i=0; i<len; i++) { pv[i]=0; sv[i]=0; zv[i]=0; }

    #pragma omp parallel
    {
        gsl_matrix *Yb = gsl_matrix_alloc(n, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
        gsl_matrix *Bb = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
        #pragma omp for schedule(dynamic)
        for(size_t ss = 0; ss < m; ss += SAMPLE_STRIP_SIZE) {
            size_t cs = (ss + SAMPLE_STRIP_SIZE > m) ? (m - ss) : SAMPLE_STRIP_SIZE;
            for(int bs = 0; bs < nrand; bs += PERM_BATCH_SIZE) {
                int cp = (bs + PERM_BATCH_SIZE > nrand) ? (nrand - bs) : PERM_BATCH_SIZE;
                for(int ip = 0; ip < cp; ip++) {
                    int *pidx = &p_tab[(size_t)(bs + ip) * n];
                    for(size_t sl = 0; sl < cs; sl++) {
                        double *src = gsl_matrix_ptr(&Yt.matrix, ss + sl, 0);
                        for(size_t g = 0; g < n; g++) gsl_matrix_set(Yb, g, (ip * cs) + sl, src[pidx[g]]);
                    }
                }
                gsl_matrix_view Ysub = gsl_matrix_submatrix(Yb, 0, 0, n, cp * cs);
                gsl_matrix_view Bsub = gsl_matrix_submatrix(Bb, 0, 0, p, cp * cs);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Ysub.matrix, 0.0, &Bsub.matrix);
                for(int ip = 0; ip < cp; ip++) {
                    for(size_t r = 0; r < p; r++) {
                        for(size_t sl = 0; sl < cs; sl++) {
                            R_xlen_t idx = (R_xlen_t)r * m + (ss + sl);
                            double brv = gsl_matrix_get(&Bsub.matrix, r, (ip * cs) + sl);
                            if(fabs(brv) >= fabs(bv[idx])) pv[idx] += 1.0;
                            zv[idx] += brv; sv[idx] += (brv * brv);
                        }
                    }
                }
            }
        }
        gsl_matrix_free(Yb); gsl_matrix_free(Bb);
    }
    double inv_n = 1.0 / nrand;
    #pragma omp parallel for schedule(static)
    for(R_xlen_t i=0; i<len; i++) {
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
        double mr = zv[i] * inv_n;
        double vr = (sv[i] * inv_n) - (mr * mr);
        double sr = sqrt(fmax(0, vr));
        sv[i] = sr; zv[i] = (sr > 1e-12) ? (bv[i] - mr) / sr : 0.0;
    }
    free(p_tab); gsl_matrix_free(I); gsl_matrix_free(T);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 3: NEW2 (Multi-threaded, T Column Permutation)
   ========================================================= */
SEXP ridgeRegTperm_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s, SEXP ncores_s) {
    int ncores = asInteger(ncores_s);
    #ifdef _OPENMP
    if (ncores > 0) omp_set_num_threads(ncores);
    #endif

    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_s, R_DimSymbol);
    size_t m = (size_t)INTEGER(y_dim)[1];
    double lambda = asReal(lambda_s);
    int nrand = asInteger(nrand_s);

    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len));
    SEXP s_s = PROTECT(allocVector(REALSXP, len));
    SEXP z_s = PROTECT(allocVector(REALSXP, len));
    SEXP p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_s), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_s), m, n);
    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I, &Xt.matrix, 0.0, T);
    
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, &Yt.matrix, 0.0, &bo.matrix);

    gsl_matrix *Tt_o = gsl_matrix_alloc(n, p);
    gsl_matrix_transpose_memcpy(Tt_o, T);

    int *p_tab = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *tmp = (int*)malloc(n * sizeof(int));
    for(size_t k=0; k<n; k++) tmp[k] = (int)k;
    srand(0);
    for(int i=0; i<nrand; i++) { shuffle(tmp, (int)n); memcpy(&p_tab[(size_t)i * n], tmp, n * sizeof(int)); }
    free(tmp);

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
        double *ms = &th_s[(size_t)tid * len], *mq = &th_q[(size_t)tid * len], *mc = &th_c[(size_t)tid * len];
        gsl_matrix *Ttp = gsl_matrix_alloc(n, p);
        gsl_matrix *br = gsl_matrix_alloc(p, m);
        size_t rsz = p * sizeof(double);
        #pragma omp for schedule(dynamic, 16)
        for(int i_r = 0; i_r < nrand; i_r++) {
            int *pidx = &p_tab[(size_t)i_r * n];
            for(size_t i = 0; i < n; i++) memcpy(Ttp->data + (pidx[i] * p), Tt_o->data + (i * p), rsz);
            gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Ttp, &Yt.matrix, 0.0, br);
            for(size_t r = 0; r < p; r++) {
                for(size_t c = 0; c < m; c++) {
                    R_xlen_t idx = (R_xlen_t)r * m + c;
                    double brv = gsl_matrix_get(br, r, c);
                    ms[idx] += brv; mq[idx] += brv * brv;
                    if(fabs(brv) >= fabs(bv[idx])) mc[idx] += 1.0;
                }
            }
        }
        gsl_matrix_free(Ttp); gsl_matrix_free(br);
    }

    #pragma omp parallel for schedule(static)
    for(R_xlen_t i = 0; i < len; i++) {
        double s = 0, sq = 0, c = 0;
        for(int t = 0; t < actual_threads; t++) {
            size_t off = (size_t)t * len + i;
            s += th_s[off]; sq += th_q[off]; c += th_c[off];
        }
        zv[i] = s; sv[i] = sq; pv[i] = c;
    }
    free(th_s); free(th_q); free(th_c); free(p_tab);
    double inv_n = 1.0 / nrand;
    #pragma omp parallel for schedule(static)
    for(R_xlen_t i = 0; i < len; i++) {
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
        double mr = zv[i] * inv_n;
        double vr = (sv[i] * inv_n) - (mr * mr);
        double sr = sqrt(fmax(0, vr));
        sv[i] = sr; zv[i] = (sr > 1e-12) ? (bv[i] - mr) / sr : 0.0;
    }
    gsl_matrix_free(Tt_o); gsl_matrix_free(I); gsl_matrix_free(T);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   REGISTRATION
   ========================================================= */
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
