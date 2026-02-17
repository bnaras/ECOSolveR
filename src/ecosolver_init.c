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

extern SEXP ecos_setup_R(SEXP sMNP, SEXP sl, SEXP sq, SEXP se,
			 SEXP sGpr, SEXP sGjc, SEXP sGir,
			 SEXP sApr, SEXP sAjc, SEXP sAir,
			 SEXP sc, SEXP sh, SEXP sb,
			 SEXP ecosControl);

extern SEXP ecos_solve_R(SEXP ext_ptr, SEXP ecosControl);

extern SEXP ecos_update_R(SEXP ext_ptr, SEXP sGpr, SEXP sApr,
			  SEXP sc, SEXP sh, SEXP sb);

extern SEXP ecos_cleanup_R(SEXP ext_ptr);

static const R_CallMethodDef CallEntries[] = {
  {"ecos_csolve",    (DL_FUNC) &ecos_csolve,    16},
  {"ecos_setup_R",   (DL_FUNC) &ecos_setup_R,   14},
  {"ecos_solve_R",   (DL_FUNC) &ecos_solve_R,    2},
  {"ecos_update_R",  (DL_FUNC) &ecos_update_R,   6},
  {"ecos_cleanup_R", (DL_FUNC) &ecos_cleanup_R,  1},
  {NULL, NULL, 0}
};

void R_init_ECOSolveR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
