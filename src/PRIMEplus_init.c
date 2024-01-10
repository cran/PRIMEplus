#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/


/* .C calls */
extern void C_ReRandLRT(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_LRT(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_EM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"C_ReRandLRT", (DL_FUNC) &C_ReRandLRT, 14},
    {"C_LRT",       (DL_FUNC) &C_LRT,       16},
    {"C_EM",        (DL_FUNC) &C_EM,        15},
    {NULL, NULL, 0}
};

void R_init_PRIMEplus(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
