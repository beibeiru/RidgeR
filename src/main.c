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

gsl_matrix *RVectorObject_to_gsl_matrix(double *vec, size_t nr, size_t nc) {
    gsl_block *b = (gsl_block *)malloc(sizeof(gsl_block));
    gsl_matrix *r = (gsl_matrix *)malloc(sizeof(gsl_matrix));
    r->size1 = nr;
    r->tda = r->size2 = nc;
    r->owner = 1;
    b->data = r->data = vec;
    r->block = b;
    b->size = r->size1 * r->size2;
    return r;
}

void gsl_matrix_partial_free(gsl_matrix *x) {
    if (x) {
        if (x->block) free(x->block);
        free(x);
    }
}

void shuffle(int array[], const int n) {
    int i, j, t;
    for (i = 0; i < n - 1; i++) {
        j = i + rand() / (RAND_MAX / (n - i) + 1);
        t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

// Internal helper to wrap results into a named R list
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

    gsl_matrix *X = RVectorObject_to_gsl_matrix(REAL(X_sexp), n, p);
    gsl_matrix *Y = RVectorObject_to_gsl_matrix(REAL(Y_sexp), n, m);
    gsl_matrix *beta = RVectorObject_to_gsl_matrix(REAL(beta_s), p, m);
    gsl_matrix *aver_sq = RVectorObject_to_gsl_matrix(REAL(se_s), p, m);
    gsl_matrix *zscore = RVectorObject_to_gsl_matrix(REAL(zs_s), p, m);
    gsl_matrix *pvalue = RVectorObject_to_gsl_matrix(REAL(pv_s), p, m);

    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y_rand = gsl_matrix_alloc(n, m);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
    gsl_matrix *aver = gsl_matrix_alloc(p, m);

    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

    int *array_index = (int*)malloc(n * sizeof(int));
    for(int i=0; i<n; i++) array_index[i] = i;
    srand(0);
    gsl_matrix_set_zero(aver); gsl_matrix_set_zero(aver_sq); gsl_matrix_set_zero(pvalue);

    for(int i=0; i<nrand; i++) {
        shuffle(array_index, n);
        for(int j=0; j<n; j++){
            gsl_vector_const_view t_v = gsl_matrix_const_row(Y, array_index[j]);
            gsl_matrix_set_row(Y_rand, j, &t_v.vector);
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);
        for(int j=0; j < p*m; j++) {
            if(fabs(beta_rand->data[j]) >= fabs(beta->data[j])) pvalue->data[j]++;
        }
        gsl_matrix_add(aver, beta_rand);
        gsl_matrix_mul_elements(beta_rand, beta_rand);
        gsl_matrix_add(aver_sq, beta_rand);
    }

    gsl_matrix_scale(aver, 1.0/nrand);
    gsl_matrix_scale(aver_sq, 1.0/nrand);
    gsl_matrix_add_constant(pvalue, 1.0);
    gsl_matrix_scale(pvalue, 1.0/(nrand + 1.0));
    gsl_matrix_memcpy(zscore, beta);
    gsl_matrix_sub(zscore, aver);
    gsl_matrix_mul_elements(aver, aver);
    gsl_matrix_sub(aver_sq, aver);
    for(int i=0; i < p*m; i++) aver_sq->data[i] = sqrt(fmax(0, aver_sq->data[i]));
    gsl_matrix_div_elements(zscore, aver_sq);

    free(array_index);
    gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_rand); gsl_matrix_free(beta_rand); gsl_matrix_free(aver);
    gsl_matrix_partial_free(X); gsl_matrix_partial_free(Y);
    gsl_matrix_partial_free(beta); gsl_matrix_partial_free(aver_sq);
    gsl_matrix_partial_free(zscore); gsl_matrix_partial_free(pvalue);

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

    gsl_matrix *X = RVectorObject_to_gsl_matrix(REAL(X_sexp), n, p);
    gsl_matrix *Y = RVectorObject_to_gsl_matrix(REAL(Y_sexp), n, m);
    gsl_matrix *beta = RVectorObject_to_gsl_matrix(REAL(beta_s), p, m);

    gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix_set_identity(I_mat);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I_mat);
    gsl_linalg_cholesky_decomp(I_mat);
    gsl_linalg_cholesky_invert(I_mat);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I_mat, X, 0.0, T);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

    gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
    gsl_matrix_transpose_memcpy(Tt_orig, T);
    gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);

    int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *temp_idx = (int*)malloc(n * sizeof(int));
    for(int k=0; k<n; k++) temp_idx[k] = k;
    srand(0);
    for(int i_rand = 0; i_rand < nrand; i_rand++) {
        shuffle(temp_idx, n);
        memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
    }
    free(temp_idx);

    double *pv_p = REAL(pv_s), *se_p = REAL(se_s), *zs_p = REAL(zs_s), *b_p = REAL(beta_s);
    for(size_t i=0; i<total_len; i++) { pv_p[i]=0; se_p[i]=0; zs_p[i]=0; }

    size_t row_size = p * sizeof(double);
    for(int i_rand=0; i_rand<nrand; i_rand++) {
        int *p_idx = &perm_table[(size_t)i_rand * n];
        for(int i=0; i<n; i++) memcpy(Tt_perm->data + (p_idx[i] * p), Tt_orig->data + (i * p), row_size);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Tt_perm, Y, 0.0, beta_rand);
        for(size_t k=0; k<total_len; k++) {
            double b_rnd = beta_rand->data[k];
            if(fabs(b_rnd) >= fabs(b_p[k])) pv_p[k] += 1.0;
            zs_p[k] += b_rnd; se_p[k] += b_rnd * b_rnd;
        }
    }

    double inv_nr = 1.0 / nrand;
    for(size_t i=0; i<total_len; i++) {
        pv_p[i] = (pv_p[i] + 1.0) / (nrand + 1.0);
        double mean_r = zs_p[i] * inv_nr;
        double var_r = (se_p[i] * inv_nr) - (mean_r * mean_r);
        double se_r = sqrt(fmax(0, var_r));
        se_p[i] = se_r;
        zs_p[i] = (se_r > 1e-12) ? (b_p[i] - mean_r) / se_r : 0.0;
    }

    free(perm_table); gsl_matrix_free(Tt_orig); gsl_matrix_free(Tt_perm); gsl_matrix_free(beta_rand);
    gsl_matrix_free(I_mat); gsl_matrix_free(T);
    gsl_matrix_partial_free(X); gsl_matrix_partial_free(Y); gsl_matrix_partial_free(beta);

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

    gsl_matrix *Xt = RVectorObject_to_gsl_matrix(REAL(X_sexp), p, n);
    gsl_matrix *Yt = RVectorObject_to_gsl_matrix(REAL(Y_sexp), m, n);
    gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I_mat);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
    gsl_linalg_cholesky_decomp(I_mat);
    gsl_linalg_cholesky_invert(I_mat);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);
    
    gsl_matrix *beta_obs = RVectorObject_to_gsl_matrix(b_vec, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta_obs);

    int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *temp_idx = (int*)malloc(n * sizeof(int));
    for(size_t k=0; k<n; k++) temp_idx[k] = (int)k;
    srand(0);
    for(int i_rand = 0; i_rand < nrand; i_rand++) {
        shuffle(temp_idx, (int)n);
        memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
    }
    free(temp_idx);

    #pragma omp parallel for schedule(static)
    for(R_xlen_t i=0; i<total_len; i++) { pv_vec[i]=0; se_vec[i]=0; zs_vec[i]=0; }

    #pragma omp parallel
    {
        gsl_matrix *Y_block = gsl_matrix_alloc(n, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);
        gsl_matrix *B_block = gsl_matrix_alloc(p, SAMPLE_STRIP_SIZE * PERM_BATCH_SIZE);

        #pragma omp for schedule(dynamic)
        for(size_t samp_start = 0; samp_start < m; samp_start += SAMPLE_STRIP_SIZE) {
            size_t curr_s = (samp_start + SAMPLE_STRIP_SIZE > m) ? (m - samp_start) : SAMPLE_STRIP_SIZE;
            for(int b_start = 0; b_start < nrand; b_start += PERM_BATCH_SIZE) {
                int curr_p = (b_start + PERM_BATCH_SIZE > nrand) ? (nrand - b_start) : PERM_BATCH_SIZE;
                for(int i_p = 0; i_p < curr_p; i_p++) {
                    int *p_idx = &perm_table[(size_t)(b_start + i_p) * n];
                    for(size_t s_l = 0; s_l < curr_s; s_l++) {
                        double *src = gsl_matrix_ptr(Yt, samp_start + s_l, 0);
                        for(size_t g = 0; g < n; g++) gsl_matrix_set(Y_block, g, (i_p * curr_s) + s_l, src[p_idx[g]]);
                    }
                }
                gsl_matrix_view Y_sub = gsl_matrix_submatrix(Y_block, 0, 0, n, curr_p * curr_s);
                gsl_matrix_view B_sub = gsl_matrix_submatrix(B_block, 0, 0, p, curr_p * curr_s);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y_sub.matrix, 0.0, &B_sub.matrix);
                for(int i_p = 0; i_p < curr_p; i_p++) {
                    for(size_t r = 0; r < p; r++) {
                        for(size_t s_l = 0; s_l < curr_s; s_l++) {
                            R_xlen_t idx = (R_xlen_t)r * m + (samp_start + s_l);
                            double b_rnd = gsl_matrix_get(&B_sub.matrix, r, (i_p * curr_s) + s_l);
                            if(fabs(b_rnd) >= fabs(b_vec[idx])) pv_vec[idx] += 1.0;
                            zs_vec[idx] += b_rnd; se_vec[idx] += (b_rnd * b_rnd);
                        }
                    }
                }
            }
        }
        gsl_matrix_free(Y_block); gsl_matrix_free(B_block);
    }

    double inv_nrand = 1.0 / nrand;
    #pragma omp parallel for schedule(static)
    for(R_xlen_t i=0; i<total_len; i++) {
        pv_vec[i] = (pv_vec[i] + 1.0) / (nrand + 1.0);
        double mean_rand = zs_vec[i] * inv_nrand;
        double var_rand = (se_vec[i] * inv_nrand) - (mean_rand * mean_rand);
        double se_rand = sqrt(fmax(0, var_rand));
        se_vec[i] = se_rand;
        zs_vec[i] = (se_rand > 1e-12) ? (b_vec[i] - mean_rand) / se_rand : 0.0;
    }

    free(perm_table); gsl_matrix_free(I_mat); gsl_matrix_free(T);
    gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta_obs);
    
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

    gsl_matrix *Xt = RVectorObject_to_gsl_matrix(REAL(X_sexp), p, n);
    gsl_matrix *Yt = RVectorObject_to_gsl_matrix(REAL(Y_sexp), m, n);
    gsl_matrix *I_mat = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I_mat);
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, Xt, lambda, I_mat);
    gsl_linalg_cholesky_decomp(I_mat);
    gsl_linalg_cholesky_invert(I_mat);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, I_mat, Xt, 0.0, T);
    
    gsl_matrix *beta_obs = RVectorObject_to_gsl_matrix(b_vec, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, T, Yt, 0.0, beta_obs);

    gsl_matrix *Tt_orig = gsl_matrix_alloc(n, p);
    gsl_matrix_transpose_memcpy(Tt_orig, T);

    int *perm_table = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *temp_idx = (int*)malloc(n * sizeof(int));
    for(size_t k=0; k<n; k++) temp_idx[k] = (int)k;
    srand(0);
    for(int i_rand = 0; i_rand < nrand; i_rand++) {
        shuffle(temp_idx, (int)n);
        memcpy(&perm_table[(size_t)i_rand * n], temp_idx, n * sizeof(int));
    }
    free(temp_idx);

    int actual_threads = 1;
    #ifdef _OPENMP
    #pragma omp parallel
    { #pragma omp single actual_threads = omp_get_num_threads(); }
    #endif

    double *th_sum = (double*)calloc((size_t)actual_threads * total_len, sizeof(double));
    double *th_sq  = (double*)calloc((size_t)actual_threads * total_len, sizeof(double));
    double *th_cnt = (double*)calloc((size_t)actual_threads * total_len, sizeof(double));

    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        double *my_s = &th_sum[(size_t)tid * total_len], *my_sq = &th_sq[(size_t)tid * total_len], *my_c = &th_cnt[(size_t)tid * total_len];
        gsl_matrix *Tt_perm = gsl_matrix_alloc(n, p);
        gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
        size_t r_bytes = p * sizeof(double);
        #pragma omp for schedule(dynamic, 16)
        for(int i_rand = 0; i_rand < nrand; i_rand++) {
            int *p_idx = &perm_table[(size_t)i_rand * n];
            for(size_t i = 0; i < n; i++) memcpy(Tt_perm->data + (p_idx[i] * p), Tt_orig->data + (i * p), r_bytes);
            gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, Tt_perm, Yt, 0.0, beta_rand);
            for(size_t r = 0; r < p; r++) {
                for(size_t c = 0; c < m; c++) {
                    R_xlen_t idx = (R_xlen_t)r * m + c;
                    double b_rnd = gsl_matrix_get(beta_rand, r, c);
                    my_s[idx] += b_rnd; my_sq[idx] += b_rnd * b_rnd;
                    if(fabs(b_rnd) >= fabs(b_vec[idx])) my_c[idx] += 1.0;
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
            s += th_sum[off]; sq += th_sq[off]; c += th_cnt[off];
        }
        zs_vec[i] = s; se_vec[i] = sq; pv_vec[i] = c;
    }

    free(th_sum); free(th_sq); free(th_cnt); free(perm_table);
    double inv_nrand = 1.0 / nrand;
    #pragma omp parallel for schedule(static)
    for(R_xlen_t i = 0; i < total_len; i++) {
        pv_vec[i] = (pv_vec[i] + 1.0) / (nrand + 1.0);
        double mean_r = zs_vec[i] * inv_nrand;
        double var_r = (se_vec[i] * inv_nrand) - (mean_r * mean_r);
        double se_r = sqrt(fmax(0, var_r));
        se_vec[i] = se_r;
        zs_vec[i] = (se_r > 1e-12) ? (b_vec[i] - mean_r) / se_r : 0.0;
    }

    gsl_matrix_free(Tt_orig); gsl_matrix_free(I_mat); gsl_matrix_free(T);
    gsl_matrix_partial_free(Xt); gsl_matrix_partial_free(Yt); gsl_matrix_partial_free(beta_obs);
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
