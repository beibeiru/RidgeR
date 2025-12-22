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

#define SAMPLE_STRIP_SIZE 64
#define PERM_BATCH_SIZE 64

/* =========================================================
   HELPER FUNCTIONS
   ========================================================= */

// Standard shuffle
void shuffle(int array[], const int n) {
    int i, j, t;
    for (i = 0; i < n - 1; i++) {
        j = i + rand() / (RAND_MAX / (n - i) + 1);
        t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

// Wrapper for R results list
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
   VERSION 1: OLD (Single-threaded, Y row permutation)
   ========================================================= */
SEXP ridgeReg_old_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp) {
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    int n = INTEGER(x_dim)[0], p = INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
    int m = INTEGER(y_dim)[1];
    int nrand = asInteger(nrand_sexp);
    double lambda = asReal(lambda_sexp);

    R_xlen_t total_len = (R_xlen_t)p * m;
    SEXP beta_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP se_s   = PROTECT(allocVector(REALSXP, total_len));
    SEXP zs_s   = PROTECT(allocVector(REALSXP, total_len));
    SEXP pv_s   = PROTECT(allocVector(REALSXP, total_len));

    gsl_matrix_view X = gsl_matrix_view_array(REAL(X_sexp), n, p);
    gsl_matrix_view Y = gsl_matrix_view_array(REAL(Y_sexp), n, m);
    gsl_matrix_view beta = gsl_matrix_view_array(REAL(beta_s), p, m);
    gsl_matrix_view aver_sq = gsl_matrix_view_array(REAL(se_s), p, m);
    gsl_matrix_view zscore = gsl_matrix_view_array(REAL(zs_s), p, m);
    gsl_matrix_view pvalue = gsl_matrix_view_array(REAL(pv_s), p, m);

    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y_rand = gsl_matrix_alloc(n, m);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
    gsl_matrix *aver = gsl_matrix_alloc(p, m);

    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &X.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, &X.matrix, 0.0, T);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y.matrix, 0.0, &beta.matrix);

    int *idx = (int*)malloc(n * sizeof(int));
    for(int i=0; i<n; i++) idx[i] = i;
    srand(0);
    gsl_matrix_set_zero(aver); gsl_matrix_set_zero(&aver_sq.matrix); gsl_matrix_set_zero(&pvalue.matrix);

    for(int i=0; i<nrand; i++) {
        shuffle(idx, n);
        for(int j=0; j<n; j++){
            gsl_vector_const_view t_v = gsl_matrix_const_row(&Y.matrix, idx[j]);
            gsl_matrix_set_row(Y_rand, j, &t_v.vector);
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);
        for(int j=0; j < p*m; j++) {
            if(fabs(beta_rand->data[j]) >= fabs(beta.matrix.data[j])) pvalue.matrix.data[j]++;
        }
        gsl_matrix_add(aver, beta_rand);
        gsl_matrix_mul_elements(beta_rand, beta_rand);
        gsl_matrix_add(&aver_sq.matrix, beta_rand);
    }

    gsl_matrix_scale(aver, 1.0/nrand);
    gsl_matrix_scale(&aver_sq.matrix, 1.0/nrand);
    gsl_matrix_add_constant(&pvalue.matrix, 1.0);
    gsl_matrix_scale(&pvalue.matrix, 1.0/(nrand + 1.0));
    gsl_matrix_memcpy(&zscore.matrix, &beta.matrix);
    gsl_matrix_sub(&zscore.matrix, aver);
    gsl_matrix_mul_elements(aver, aver);
    gsl_matrix_sub(&aver_sq.matrix, aver);
    for(int i=0; i < p*m; i++) aver_sq.matrix.data[i] = sqrt(fmax(0, aver_sq.matrix.data[i]));
    gsl_matrix_div_elements(&zscore.matrix, &aver_sq.matrix);

    free(idx);
    gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_rand); gsl_matrix_free(beta_rand); gsl_matrix_free(aver);

    SEXP res = create_res_list(beta_s, se_s, zs_s, pv_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 1.5: OLD2 (Single-threaded, T column permutation)
   ========================================================= */
SEXP ridgeRegTperm_old_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp) {
    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    int n = INTEGER(x_dim)[0], p = INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
    int m = INTEGER(y_dim)[1];
    int nrand = asInteger(nrand_sexp);
    double lambda = asReal(lambda_sexp);

    R_xlen_t total_len = (R_xlen_t)p * m;
    SEXP beta_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP se_s   = PROTECT(allocVector(REALSXP, total_len));
    SEXP zs_s   = PROTECT(allocVector(REALSXP, total_len));
    SEXP pv_s   = PROTECT(allocVector(REALSXP, total_len));

    gsl_matrix_view X = gsl_matrix_view_array(REAL(X_sexp), n, p);
    gsl_matrix_view Y = gsl_matrix_view_array(REAL(Y_sexp), n, m);
    gsl_matrix_view beta = gsl_matrix_view_array(REAL(beta_s), p, m);

    gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix_set_identity(I_mat);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &X.matrix, lambda, I_mat);
    gsl_linalg_cholesky_decomp(I_mat);
    gsl_linalg_cholesky_invert(I_mat);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I_mat, &X.matrix, 0.0, T);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y.matrix, 0.0, &beta.matrix);

    gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
    gsl_matrix_transpose_memcpy(Tt_orig, T);
    gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

    int *p_tab = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *tmp = (int*)malloc(n * sizeof(int));
    for(int k=0; k<n; k++) tmp[k] = k;
    srand(0);
    for(int i=0; i<nrand; i++) {
        shuffle(tmp, n);
        memcpy(&p_tab[(size_t)i * n], tmp, n * sizeof(int));
    }
    free(tmp);

    double *pv_p = REAL(pv_s), *se_p = REAL(se_s), *zs_p = REAL(zs_s), *b_p = REAL(beta_s);
    for(size_t i=0; i<total_len; i++) { pv_p[i]=0; se_p[i]=0; zs_p[i]=0; }

    size_t row_sz = p * sizeof(double);
    for(int i_r=0; i_r<nrand; i_r++) {
        int *p_idx = &p_tab[(size_t)i_r * n];
        for(int i=0; i<n; i++) memcpy(Tt_perm->data + (p_idx[i] * p), Tt_orig->data + (i * p), row_sz);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Tt_perm, &Y.matrix, 0.0, beta_rand);
        for(size_t k=0; k<total_len; k++) {
            double br = beta_rand->data[k];
            if(fabs(br) >= fabs(b_p[k])) pv_p[k] += 1.0;
            zs_p[k] += br; se_p[k] += br * br;
        }
    }

    double inv_n = 1.0 / nrand;
    for(size_t i=0; i<total_len; i++) {
        pv_p[i] = (pv_p[i] + 1.0) / (nrand + 1.0);
        double mr = zs_p[i] * inv_n;
        double vr = (se_p[i] * inv_n) - (mr * mr);
        double sr = sqrt(fmax(0, vr));
        se_p[i] = sr;
        zs_p[i] = (sr > 1e-12) ? (b_p[i] - mr) / sr : 0.0;
    }

    free(p_tab); gsl_matrix_free(Tt_orig); gsl_matrix_free(Tt_perm); gsl_matrix_free(beta_rand);
    gsl_matrix_free(I_mat); gsl_matrix_free(T);

    SEXP res = create_res_list(beta_s, se_s, zs_s, pv_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 2: NEW (Multi-threaded, Y row permutation)
   ========================================================= */
SEXP ridgeRegFast_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp, SEXP ncores_sexp) {
    int ncores = asInteger(ncores_sexp);
    #ifdef _OPENMP
    if (ncores > 0) omp_set_num_threads(ncores);
    #endif

    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
    size_t m = (size_t)INTEGER(y_dim)[1];
    double lambda = asReal(lambda_sexp);
    int nrand = asInteger(nrand_sexp);

    R_xlen_t total_len = (R_xlen_t)p * m;
    SEXP beta_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP se_s   = PROTECT(allocVector(REALSXP, total_len));
    SEXP zs_s   = PROTECT(allocVector(REALSXP, total_len));
    SEXP pv_s   = PROTECT(allocVector(REALSXP, total_len));
    double *b_vec = REAL(beta_s), *se_vec = REAL(se_s), *zs_vec = REAL(zs_s), *pv_vec = REAL(pv_s);

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);
    gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I_mat);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I_mat);
    gsl_linalg_cholesky_decomp(I_mat);
    gsl_linalg_cholesky_invert(I_mat);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, &Xt.matrix, 0.0, T);
    
    gsl_matrix_view b_obs_v = gsl_matrix_view_array(b_vec, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, &Yt.matrix, 0.0, &b_obs_v.matrix);

    int *p_tab = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *tmp = (int*)malloc(n * sizeof(int));
    for(size_t k=0; k<n; k++) tmp[k] = (int)k;
    srand(0);
    for(int i=0; i<nrand; i++) {
        shuffle(tmp, (int)n);
        memcpy(&p_tab[(size_t)i * n], tmp, n * sizeof(int));
    }
    free(tmp);

    #pragma omp parallel for schedule(static)
    for(R_xlen_t i=0; i<total_len; i++) { pv_vec[i]=0; se_vec[i]=0; zs_vec[i]=0; }

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
                gsl_matrix_view Y_sub = gsl_matrix_submatrix(Yb, 0, 0, n, cp * cs);
                gsl_matrix_view B_sub = gsl_matrix_submatrix(Bb, 0, 0, p, cp * cs);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);
                for(int ip = 0; ip < cp; ip++) {
                    for(size_t r = 0; r < p; r++) {
                        for(size_t sl = 0; sl < cs; sl++) {
                            R_xlen_t idx = (R_xlen_t)r * m + (ss + sl);
                            double br = gsl_matrix_get(&B_sub.matrix, r, (ip * cs) + sl);
                            if(fabs(br) >= fabs(b_vec[idx])) pv_vec[idx] += 1.0;
                            zs_vec[idx] += br; se_vec[idx] += (br * br);
                        }
                    }
                }
            }
        }
        gsl_matrix_free(Yb); gsl_matrix_free(Bb);
    }

    double inv_n = 1.0 / nrand;
    #pragma omp parallel for schedule(static)
    for(R_xlen_t i=0; i<total_len; i++) {
        pv_vec[i] = (pv_vec[i] + 1.0) / (nrand + 1.0);
        double mr = zs_vec[i] * inv_n;
        double vr = (se_vec[i] * inv_n) - (mr * mr);
        double sr = sqrt(fmax(0, vr));
        se_vec[i] = sr;
        zs_vec[i] = (sr > 1e-12) ? (b_vec[i] - mr) / sr : 0.0;
    }

    free(p_tab); gsl_matrix_free(I_mat); gsl_matrix_free(T);
    
    SEXP res = create_res_list(beta_s, se_s, zs_s, pv_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 3: NEW2 (Multi-threaded, T column permutation)
   ========================================================= */
SEXP ridgeRegTperm_interface(SEXP X_sexp, SEXP Y_sexp, SEXP lambda_sexp, SEXP nrand_sexp, SEXP ncores_sexp) {
    int ncores = asInteger(ncores_sexp);
    #ifdef _OPENMP
    if (ncores > 0) omp_set_num_threads(ncores);
    #endif

    SEXP x_dim = getAttrib(X_sexp, R_DimSymbol);
    size_t n = (size_t)INTEGER(x_dim)[0], p = (size_t)INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_sexp, R_DimSymbol);
    size_t m = (size_t)INTEGER(y_dim)[1];
    double lambda = asReal(lambda_sexp);
    int nrand = asInteger(nrand_sexp);

    R_xlen_t total_len = (R_xlen_t)p * m;
    SEXP beta_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP se_s   = PROTECT(allocVector(REALSXP, total_len));
    SEXP zs_s   = PROTECT(allocVector(REALSXP, total_len));
    SEXP pv_s   = PROTECT(allocVector(REALSXP, total_len));
    double *b_vec = REAL(beta_s), *se_vec = REAL(se_s), *zs_vec = REAL(zs_s), *pv_vec = REAL(pv_s);

    gsl_matrix_view Xt = gsl_matrix_view_array(REAL(X_sexp), p, n);
    gsl_matrix_view Yt = gsl_matrix_view_array(REAL(Y_sexp), m, n);
    gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I_mat);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &Xt.matrix, lambda, I_mat);
    gsl_linalg_cholesky_decomp(I_mat);
    gsl_linalg_cholesky_invert(I_mat);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, &Xt.matrix, 0.0, T);
    
    gsl_matrix_view b_obs_v = gsl_matrix_view_array(b_vec, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, &Yt.matrix, 0.0, &b_obs_v.matrix);

    gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
    gsl_matrix_transpose_memcpy(Tt_orig, T);

    int *p_tab = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *tmp = (int*)malloc(n * sizeof(int));
    for(size_t k=0; k<n; k++) tmp[k] = (int)k;
    srand(0);
    for(int i=0; i<nrand; i++) {
        shuffle(tmp, (int)n);
        memcpy(&p_tab[(size_t)i * n], tmp, n * sizeof(int));
    }
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

    double *th_s = (double*)calloc((size_t)actual_threads * total_len, sizeof(double));
    double *th_q = (double*)calloc((size_t)actual_threads * total_len, sizeof(double));
    double *th_c = (double*)calloc((size_t)actual_threads * total_len, sizeof(double));

    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        double *my_s = &th_s[(size_t)tid * total_len], *my_q = &th_q[(size_t)tid * total_len], *my_c = &th_c[(size_t)tid * total_len];
        gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);
        gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
        size_t r_sz = p * sizeof(double);
        #pragma omp for schedule(dynamic, 16)
        for(int i_r = 0; i_r < nrand; i_r++) {
            int *pidx = &p_tab[(size_t)i_r * n];
            for(size_t i = 0; i < n; i++) memcpy(Tt_perm->data + (pidx[i] * p), Tt_orig->data + (i * p), r_sz);
            gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Tt_perm, &Yt.matrix, 0.0, beta_rand);
            for(size_t r = 0; r < p; r++) {
                for(size_t c = 0; c < m; c++) {
                    R_xlen_t idx = (R_xlen_t)r * m + c;
                    double br = gsl_matrix_get(beta_rand, r, c);
                    my_s[idx] += br; my_q[idx] += br * br;
                    if(fabs(br) >= fabs(b_vec[idx])) my_c[idx] += 1.0;
                }
            }
        }
        gsl_matrix_free(Tt_perm); gsl_matrix_free(beta_rand);
    }

    #pragma omp parallel for schedule(static)
    for(R_xlen_t i = 0; i < total_len; i++) {
        double s = 0, sq = 0, c = 0;
        for(int t = 0; t < actual_threads; t++) {
            size_t off = (size_t)t * total_len + i;
            s += th_s[off]; sq += th_q[off]; c += th_c[off];
        }
        zs_vec[i] = s; se_vec[i] = sq; pv_vec[i] = c;
    }

    free(th_s); free(th_q); free(th_c); free(p_tab);
    double inv_n = 1.0 / nrand;
    #pragma omp parallel for schedule(static)
    for(R_xlen_t i = 0; i < total_len; i++) {
        pv_vec[i] = (pv_vec[i] + 1.0) / (nrand + 1.0);
        double mr = zs_vec[i] * inv_n;
        double vr = (se_vec[i] * inv_n) - (mr * mr);
        double sr = sqrt(fmax(0, vr));
        se_vec[i] = sr;
        zs_vec[i] = (sr > 1e-12) ? (b_vec[i] - mr) / sr : 0.0;
    }

    gsl_matrix_free(Tt_orig); gsl_matrix_free(I_mat); gsl_matrix_free(T);
    SEXP res = create_res_list(beta_s, se_s, zs_s, pv_s);
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
