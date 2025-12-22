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

gsl_matrix *RVectorObject_to_gsl_matrix(double *vec, size_t nr, size_t nc) {
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

void gsl_matrix_partial_free(gsl_matrix *x) {
    if(x) {
        if(x->block) free(x->block);
        free(x);
    }
}

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
   LEGACY .C VERSION
   ========================================================= */
void ridgeReg(
    double *X_vec, double *Y_vec, int *n_pt, int *p_pt, int *m_pt,
    double *lambda_pt, double *nrand_pt,
    double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec
) {
    int n = *n_pt, p = *p_pt, m = *m_pt, nrand = (int)*nrand_pt;
    double lambda = *lambda_pt;

    gsl_matrix *X = RVectorObject_to_gsl_matrix(X_vec, (size_t)n, (size_t)p);
    gsl_matrix *Y = RVectorObject_to_gsl_matrix(Y_vec, (size_t)n, (size_t)m);
    gsl_matrix *beta = RVectorObject_to_gsl_matrix(beta_vec, (size_t)p, (size_t)m);
    gsl_matrix *aver_sq = RVectorObject_to_gsl_matrix(se_vec, (size_t)p, (size_t)m);
    gsl_matrix *zscore = RVectorObject_to_gsl_matrix(zscore_vec, (size_t)p, (size_t)m);
    gsl_matrix *pvalue = RVectorObject_to_gsl_matrix(pvalue_vec, (size_t)p, (size_t)m);

    gsl_matrix *I = gsl_matrix_alloc(p, p);
    gsl_matrix *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y_rand = gsl_matrix_alloc(n, m);
    gsl_matrix *beta_rand = gsl_matrix_alloc(p, m);
    gsl_matrix *aver = gsl_matrix_alloc(p, m);
    int *array_index = (int*)malloc(n * sizeof(int));

    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, X, 0.0, T);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, beta);

    srand(0);
    for(int i=0; i<n; i++) array_index[i] = i;
    gsl_matrix_set_zero(aver); gsl_matrix_set_zero(aver_sq); gsl_matrix_set_zero(pvalue);

    for(int i=0; i<nrand; i++) {
        shuffle(array_index, n);
        for(int j=0; j<n; j++) {
            gsl_vector_const_view t = gsl_matrix_const_row(Y, array_index[j]);
            gsl_matrix_set_row(Y_rand, j, &t.vector);
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);
        for(int j=0; j < p*m; j++) if(fabs(beta_rand->data[j]) >= fabs(beta->data[j])) pvalue->data[j]++;
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
    for(int i=0; i < p*m; i++) aver_sq->data[i] = sqrt(fmax(0.0, aver_sq->data[i]));
    gsl_matrix_div_elements(zscore, aver_sq);

    free(array_index); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_rand);
    gsl_matrix_free(beta_rand); gsl_matrix_free(aver);
    gsl_matrix_partial_free(X); gsl_matrix_partial_free(Y); gsl_matrix_partial_free(beta);
    gsl_matrix_partial_free(aver_sq); gsl_matrix_partial_free(zscore); gsl_matrix_partial_free(pvalue);
}

/* =========================================================
   VERSION 1: OLD (.Call - Y Row Permutation)
   ========================================================= */
SEXP ridgeReg_old_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s) {
    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    int n = INTEGER(x_dim)[0], p = INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_s, R_DimSymbol);
    int m = INTEGER(y_dim)[1];
    int nrand = asInteger(nrand_s);
    double lambda = asReal(lambda_s);

    R_xlen_t total_len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, total_len)), s_s = PROTECT(allocVector(REALSXP, total_len));
    SEXP z_s = PROTECT(allocVector(REALSXP, total_len)), p_s = PROTECT(allocVector(REALSXP, total_len));

    gsl_matrix_view X = gsl_matrix_view_array(REAL(X_s), n, p);
    gsl_matrix_view Y = gsl_matrix_view_array(REAL(Y_s), n, m);
    gsl_matrix_view beta = gsl_matrix_view_array(REAL(b_s), p, m);
    gsl_matrix_view aver_sq = gsl_matrix_view_array(REAL(s_s), p, m);
    gsl_matrix_view zscore = gsl_matrix_view_array(REAL(z_s), p, m);
    gsl_matrix_view pvalue = gsl_matrix_view_array(REAL(p_s), p, m);

    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n);
    gsl_matrix *Y_rand = gsl_matrix_alloc(n, m), *beta_rand = gsl_matrix_alloc(p, m), *aver = gsl_matrix_alloc(p, m);

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
        for(int j=0; j<n; j++) {
            gsl_vector_const_view t = gsl_matrix_const_row(&Y.matrix, idx[j]);
            gsl_matrix_set_row(Y_rand, j, &t.vector);
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y_rand, 0.0, beta_rand);
        for(int j=0; j < p*m; j++) if(fabs(beta_rand->data[j]) >= fabs(beta.matrix.data[j])) pvalue.matrix.data[j]++;
        gsl_matrix_add(aver, beta_rand);
        gsl_matrix_mul_elements(beta_rand, beta_rand);
        gsl_matrix_add(&aver_sq.matrix, beta_rand);
    }

    gsl_matrix_scale(aver, 1.0/nrand); gsl_matrix_scale(&aver_sq.matrix, 1.0/nrand);
    gsl_matrix_add_constant(&pvalue.matrix, 1.0); gsl_matrix_scale(&pvalue.matrix, 1.0/(nrand + 1.0));
    gsl_matrix_memcpy(&zscore.matrix, &beta.matrix); gsl_matrix_sub(&zscore.matrix, aver);
    gsl_matrix_mul_elements(aver, aver); gsl_matrix_sub(&aver_sq.matrix, aver);
    for(int i=0; i < p*m; i++) aver_sq.matrix.data[i] = sqrt(fmax(0.0, aver_sq.matrix.data[i]));
    gsl_matrix_div_elements(&zscore.matrix, &aver_sq.matrix);

    free(idx); gsl_matrix_free(I); gsl_matrix_free(T); gsl_matrix_free(Y_rand);
    gsl_matrix_free(beta_rand); gsl_matrix_free(aver);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 1.5: OLD2 (.Call - T Column Permutation)
   ========================================================= */
SEXP ridgeRegTperm_old_interface(SEXP X_s, SEXP Y_s, SEXP lambda_s, SEXP nrand_s) {
    SEXP x_dim = getAttrib(X_s, R_DimSymbol);
    int n = INTEGER(x_dim)[0], p = INTEGER(x_dim)[1];
    SEXP y_dim = getAttrib(Y_s, R_DimSymbol);
    int m = INTEGER(y_dim)[1];
    int nrand = asInteger(nrand_s);
    double lambda = asReal(lambda_s);

    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len)), s_s = PROTECT(allocVector(REALSXP, len));
    SEXP z_s = PROTECT(allocVector(REALSXP, len)), p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);

    gsl_matrix_view X = gsl_matrix_view_array(REAL(X_s), n, p);
    gsl_matrix_view Y = gsl_matrix_view_array(REAL(Y_s), n, m);
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n);
    
    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &X.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, &X.matrix, 0.0, T);
    
    gsl_matrix_view beta_obs = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y.matrix, 0.0, &beta_obs.matrix);

    gsl_matrix *Tt_o = gsl_matrix_alloc(n, p); gsl_matrix_transpose_memcpy(Tt_o, T);
    gsl_matrix *Tt_p = gsl_matrix_alloc(n, p); gsl_matrix *br = gsl_matrix_alloc(p, m);

    int *p_tab = (int*)malloc((size_t)nrand * n * sizeof(int));
    int *tmp = (int*)malloc(n * sizeof(int));
    for(int k=0; k<n; k++) tmp[k] = k;
    srand(0);
    for(int i=0; i<nrand; i++) { shuffle(tmp, n); memcpy(&p_tab[(size_t)i * n], tmp, n * sizeof(int)); }
    free(tmp);

    for(R_xlen_t i=0; i<len; i++) { pv[i]=0; sv[i]=0; zv[i]=0; }
    size_t row_sz = p * sizeof(double);

    for(int i_r=0; i_r<nrand; i_r++) {
        int *p_idx = &p_tab[(size_t)i_r * n];
        for(int i=0; i<n; i++) memcpy(Tt_p->data + (p_idx[i] * p), Tt_o->data + (i * p), row_sz);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Tt_p, &Y.matrix, 0.0, br);
        for(R_xlen_t k=0; k<len; k++) {
            if(fabs(br->data[k]) >= fabs(bv[k])) pv[k] += 1.0;
            zv[k] += br->data[k]; sv[k] += br->data[k] * br->data[k];
        }
    }

    double inv_n = 1.0 / nrand;
    for(R_xlen_t i=0; i<len; i++) {
        pv[i] = (pv[i] + 1.0) / (nrand + 1.0);
        double mr = zv[i] * inv_n;
        double vr = (sv[i] * inv_n) - (mr * mr);
        double sr = sqrt(fmax(0.0, vr));
        sv[i] = sr; zv[i] = (sr > 1e-12) ? (bv[i] - mr) / sr : 0.0;
    }

    free(p_tab); gsl_matrix_free(Tt_o); gsl_matrix_free(Tt_p); gsl_matrix_free(br);
    gsl_matrix_free(I); gsl_matrix_free(T);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 2: NEW (Multi-threaded - Y Row Permutation)
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
    double lambda = asReal(lambda_s); int nrand = asInteger(nrand_s);

    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len)), s_s = PROTECT(allocVector(REALSXP, len));
    SEXP z_s = PROTECT(allocVector(REALSXP, len)), p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);

    gsl_matrix_view X = gsl_matrix_view_array(REAL(X_s), n, p);
    gsl_matrix_view Y = gsl_matrix_view_array(REAL(Y_s), n, m);
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &X.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, &X.matrix, 0.0, T);
    
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y.matrix, 0.0, &bo.matrix);

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
                        for(size_t gene = 0; gene < n; gene++) {
                            double val = gsl_matrix_get(&Y.matrix, pidx[gene], ss + sl);
                            gsl_matrix_set(Yb, gene, (ip * cs) + sl, val);
                        }
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
        double sr = sqrt(fmax(0.0, vr));
        sv[i] = sr; zv[i] = (sr > 1e-12) ? (bv[i] - mr) / sr : 0.0;
    }
    free(p_tab); gsl_matrix_free(I); gsl_matrix_free(T);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   VERSION 3: NEW2 (Multi-threaded - T Column Permutation)
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
    double lambda = asReal(lambda_s); int nrand = asInteger(nrand_s);

    R_xlen_t len = (R_xlen_t)p * m;
    SEXP b_s = PROTECT(allocVector(REALSXP, len)), s_s = PROTECT(allocVector(REALSXP, len));
    SEXP z_s = PROTECT(allocVector(REALSXP, len)), p_s = PROTECT(allocVector(REALSXP, len));
    double *bv = REAL(b_s), *sv = REAL(s_s), *zv = REAL(z_s), *pv = REAL(p_s);

    gsl_matrix_view X = gsl_matrix_view_array(REAL(X_s), n, p);
    gsl_matrix_view Y = gsl_matrix_view_array(REAL(Y_s), n, m);
    gsl_matrix *I = gsl_matrix_alloc(p, p), *T = gsl_matrix_alloc(p, n);

    gsl_matrix_set_identity(I);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &X.matrix, lambda, I);
    gsl_linalg_cholesky_decomp(I);
    gsl_linalg_cholesky_invert(I);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, I, &X.matrix, 0.0, T);
    
    gsl_matrix_view bo = gsl_matrix_view_array(bv, p, m);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, &Y.matrix, 0.0, &bo.matrix);

    gsl_matrix *Tt_o = gsl_matrix_alloc(n, p); gsl_matrix_transpose_memcpy(Tt_o, T);

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
        gsl_matrix *Ttp = gsl_matrix_alloc(n, p); gsl_matrix *br = gsl_matrix_alloc(p, m);
        size_t rsz = p * sizeof(double);
        #pragma omp for schedule(dynamic, 16)
        for(int i_r = 0; i_r < nrand; i_r++) {
            int *pidx = &p_tab[(size_t)i_r * n];
            for(size_t i = 0; i < n; i++) memcpy(Ttp->data + (pidx[i] * p), Tt_o->data + (i * p), rsz);
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Ttp, &Y.matrix, 0.0, br);
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
        double sr = sqrt(fmax(0.0, vr));
        sv[i] = sr; zv[i] = (sr > 1e-12) ? (bv[i] - mr) / sr : 0.0;
    }
    gsl_matrix_free(Tt_o); gsl_matrix_free(I); gsl_matrix_free(T);
    SEXP res = create_res_list(b_s, s_s, z_s, p_s);
    UNPROTECT(4); return res;
}

/* =========================================================
   REGISTRATION
   ========================================================= */
static const R_CMethodDef cMethods[] = {
    {"ridgeReg", (DL_FUNC) &ridgeReg, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef callMethods[] = {
    {"ridgeReg_old_interface",      (DL_FUNC) &ridgeReg_old_interface,      4},
    {"ridgeRegTperm_old_interface", (DL_FUNC) &ridgeRegTperm_old_interface, 4},
    {"ridgeRegFast_interface",      (DL_FUNC) &ridgeRegFast_interface,      5},
    {"ridgeRegTperm_interface",     (DL_FUNC) &ridgeRegTperm_interface,     5},
    {NULL, NULL, 0}
};

void R_init_RidgeR(DllInfo *dll) {
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
