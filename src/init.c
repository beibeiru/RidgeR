#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void ridgeReg(double *, double *, int *, int *, int *, double *, double *, double *, double *, double *, double *);

/* .Call calls */
extern SEXP ridgeRegFast_call_wrapper(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"ridgeReg", (DL_FUNC) &ridgeReg, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"ridgeRegFast_call_wrapper", (DL_FUNC) &ridgeRegFast_call_wrapper, 5},
    {NULL, NULL, 0}
};

void R_init_RidgeR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
