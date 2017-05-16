#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/* .Call entry points */

extern SEXP ecos_csolve(SEXP sMNP, SEXP sl, SEXP sq, SEXP se,
			SEXP sGpr, SEXP sGjc, SEXP sGir,
			SEXP sApr, SEXP sAjc, SEXP sAir,
			SEXP sc, SEXP sh, SEXP sb,
			SEXP sbool_vars, SEXP sint_vars,
			SEXP ecosControl);

static const R_CallMethodDef CallEntries[] = {
  {"ecos_csolve", (DL_FUNC) &ecos_csolve, 16},
  {NULL, NULL, 0}
};

void R_init_EcoSolveR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
