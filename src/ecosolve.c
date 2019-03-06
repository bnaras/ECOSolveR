#include <R.h>
#include <Rinternals.h>

#define IGNORE_UNUSED_FUNCTIONS

#include "ecos.h"
#include "ecos_bb.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("ECOSolveR", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

/** Modified from R-exts example for getting a named element */
SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

/** Set up ECOS_BB parameters */
void setup_ECOS_BB_params (settings_bb *settings, SEXP ecosControl) {
  settings->abs_tol_gap = (pfloat) (REAL(getListElement(ecosControl, "MI_ABS_EPS"))[0]);
  settings->rel_tol_gap = (pfloat) (REAL(getListElement(ecosControl, "MI_REL_EPS"))[0]);
  settings->integer_tol = (pfloat) (REAL(getListElement(ecosControl, "MI_INT_TOL"))[0]);
  settings->maxit = (idxint) (INTEGER(getListElement(ecosControl, "MI_MAX_ITERS"))[0]);
  settings->verbose = (idxint) (INTEGER(getListElement(ecosControl, "VERBOSE"))[0]);
}

/** Set up ECOS parameters */
void setup_ECOS_params (pwork *mywork, SEXP ecosControl) {
  mywork->stgs->feastol = (pfloat) (REAL(getListElement(ecosControl, "FEASTOL"))[0]);
  mywork->stgs->reltol = (pfloat) (REAL(getListElement(ecosControl, "RELTOL"))[0]);
  mywork->stgs->abstol = (pfloat) (REAL(getListElement(ecosControl, "ABSTOL"))[0]);
  mywork->stgs->feastol_inacc = (pfloat) (REAL(getListElement(ecosControl, "FEASTOL_INACC"))[0]);
  mywork->stgs->abstol_inacc = (pfloat) (REAL(getListElement(ecosControl, "ABSTOL_INACC"))[0]);
  mywork->stgs->reltol_inacc = (pfloat) (REAL(getListElement(ecosControl, "RELTOL_INACC"))[0]);
  mywork->stgs->maxit = (idxint) (INTEGER(getListElement(ecosControl, "MAXIT"))[0]);
  mywork->stgs->verbose = (idxint) (INTEGER(getListElement(ecosControl, "VERBOSE"))[0]);
}

/**
 * This code is mightily dependent on the changes made to idxint
 * and pfloat types. In R, the Matrix package uses integer vectors
 * to point to the indices and therefore they are ints, not
 * long vectors which can be over 2^32 in length. Also, for now
 * pfloat seems to be the same as double, so we are ok. Otherwise,
 * we have to fix things. */

/* Make a copy of int vector to an idxint vector */
idxint * int2idxint(SEXP intsxp) {
  int n, *source;
  idxint *result;
  n = length(intsxp);
  source = INTEGER(intsxp);
  result = (idxint *) malloc(n * sizeof(idxint));
  for (int i = 0; i < n; i++) {
    result[i] = source[i];
  }
  return(result);
}

/* Allocation pointers */
typedef struct alloc_ptr {
  idxint *q;
  idxint *Gjc;
  idxint *Gir;
  idxint *Ajc;
  idxint *Air;
  idxint *bool_vars;
  idxint *int_vars;
} alloc_ptr;

void free_allocated(alloc_ptr *allocated) {
  if (allocated->q != NULL) free(allocated->q);
  if (allocated->Gjc != NULL) free(allocated->Gjc);
  if (allocated->Gir != NULL) free(allocated->Gir);
  if (allocated->Ajc != NULL) free(allocated->Ajc);
  if (allocated->Air != NULL) free(allocated->Air);
  if (allocated->bool_vars != NULL) free(allocated->bool_vars);
  if (allocated->int_vars != NULL) free(allocated->int_vars);
}

SEXP ecos_csolve(SEXP sMNP, SEXP sl, SEXP sq, SEXP se,
		 SEXP sGpr, SEXP sGjc, SEXP sGir,
		 SEXP sApr, SEXP sAjc, SEXP sAir,
		 SEXP sc, SEXP sh, SEXP sb,
		 SEXP sbool_vars, SEXP sint_vars,
		 SEXP ecosControl) {
  // ECOS data structures
  idxint exitcode;
  idxint n, m, p, l, ncones, *q, e, *Gjc, *Gir, *Ajc, *Air, *bool_vars, *int_vars;
  pfloat *Gpr, *Apr, *c, *h, *b;
  alloc_ptr allocated;
  const char* infostring;
  int i, nProtected, numerr, mi_iterations = -1, num_bool, num_int;

  /* Since R long_vector structure is different, we have to *
   * make copies of indices into allocated structures and   *
   * free them back when we are done. This structure keeps  *
   * track of what is allocated.                            */
  allocated.q = NULL;
  allocated.Gjc = NULL;
  allocated.Gir = NULL;
  allocated.Ajc = NULL;
  allocated.Air = NULL;
  allocated.bool_vars = NULL;
  allocated.int_vars = NULL;

  n = (idxint) (INTEGER(sMNP)[1]);
  m = (idxint) (INTEGER(sMNP)[0]);
  p = (idxint) (INTEGER(sMNP)[2]);

  l = (idxint) (INTEGER(sl)[0]);
  ncones = (idxint) length(sq);
  if (Rf_isNull(sq)) {
    q = NULL;
  } else {
    q = int2idxint(sq);
    allocated.q = q;
  }

  e = (idxint) (INTEGER(se)[0]);

  Gpr = Rf_isNull(sGpr)? NULL : (pfloat *) REAL(sGpr);
  if (Rf_isNull(sGjc)) {
    Gjc = NULL;
  } else {
    Gjc = int2idxint(sGjc);
    allocated.Gjc = Gjc;
  }
  if (Rf_isNull(sGir)) {
    Gir = NULL;
  } else {
    Gir = int2idxint(sGir);
    allocated.Gir = Gir;
  }

  Apr = Rf_isNull(sApr)? NULL : (pfloat *) REAL(sApr);
  if (Rf_isNull(sAjc)) {
    Ajc = NULL;
  } else {
    Ajc = int2idxint(sAjc);
    allocated.Ajc = Ajc;
  }
  if (Rf_isNull(sAir)) {
    Air = NULL;
  } else {
    Air = int2idxint(sAir);
    allocated.Air = Air;
  }

  c = Rf_isNull(sc)? NULL : (pfloat *) REAL(sc);
  h = Rf_isNull(sh)? NULL : (pfloat *) REAL(sh);
  b = Rf_isNull(sb)? NULL : (pfloat *) REAL(sb);

  if (Rf_isNull(sbool_vars)) {
    bool_vars = NULL;
    num_bool = 0;
  } else {
    num_bool = length(sbool_vars);
    bool_vars = int2idxint(sbool_vars);
    allocated.bool_vars = bool_vars;
  }

  if (Rf_isNull(sint_vars)) {
    int_vars = NULL;
    num_int = 0;
  } else {
    num_int = length(sint_vars);
    int_vars = int2idxint(sint_vars);
    allocated.int_vars = int_vars;
  }

  ecos_bb_pwork* myecos_bb_work = NULL;
  pwork* mywork = NULL;

  if (num_bool > 0 || num_int > 0) {
    /* This calls ECOS setup function. */
    settings_bb opts_ecos_bb;
    /* Set up parameters for ECOS_BB */
    setup_ECOS_BB_params(&opts_ecos_bb, ecosControl);

    myecos_bb_work = ECOS_BB_setup(n, m, p, l, ncones, q, e, Gpr, Gjc, Gir,
				   Apr, Ajc, Air, c, h, b,
				   num_bool, bool_vars, num_int, int_vars, &opts_ecos_bb);
    if ( myecos_bb_work == NULL ) {
      /*  ERROR NEED TO FIX UP  */
      /* "Internal problem occurred in ECOS_BB while setting up the problem.\nPlease send a bug report with data to Alexander Domahidi.\nEmail: domahidi@control.ee.ethz.ch" */
      /* Free whatever was allocated, including the structure */
      free_allocated(&allocated);
      return R_NilValue;
    }

    mywork = myecos_bb_work->ecos_prob;

    /* Set up parameters for ECOS */
    setup_ECOS_params(mywork, ecosControl);

    /* Solve! */
    exitcode = ECOS_BB_solve(myecos_bb_work);
    mi_iterations = (int) myecos_bb_work->iter;

  } else {
    /* This calls ECOS setup function */
    mywork = ECOS_setup(n, m, p, l, ncones, q, e, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
    if ( mywork == NULL ) {
      /*  ERROR NEED TO FIX UP  */
      /* "Internal problem occurred in ECOS while setting up the problem.\nPlease send a bug report with data to Alexander Domahidi.\nEmail: domahidi@control.ee.ethz.ch" */
      /* Free whatever was allocated, including the structure */
      free_allocated(&allocated);
      return R_NilValue;
    }

    /* Set up parameters for ECOS */
    setup_ECOS_params(mywork, ecosControl);
    exitcode = ECOS_solve(mywork);
  }

  nProtected = 0; /* Number of SEXPs protected */

  /* create output (all data is *deep copied*) */
  /* x */
  SEXP x;
  PROTECT(x = allocVector(REALSXP, n));
  nProtected++;
  memcpy(REAL(x), mywork->x, n * sizeof(double));

  /* y */
  SEXP y;
  PROTECT(y = allocVector(REALSXP, p));
  nProtected++;
  memcpy(REAL(y), mywork->y, p * sizeof(double));

  /* s */
  SEXP s;
  PROTECT(s = allocVector(REALSXP, m));
  nProtected++;
  memcpy(REAL(s), mywork->s, m * sizeof(double));

  /* z */
  SEXP z;
  PROTECT(z = allocVector(REALSXP, m));
  nProtected++;
  memcpy(REAL(z), mywork->z, m * sizeof(double));

  /* Infostring based on return code */
  if (num_bool > 0 || num_int > 0){
    /* info dict */
    /* infostring */
    switch( exitcode ){
        case MI_OPTIMAL_SOLN:
            infostring = "Optimal branch and bound solution found";
            break;
        case MI_MAXITER_FEASIBLE_SOLN:
            infostring = "Maximum iterations reached with feasible solution found";
            break;
        case MI_MAXITER_NO_SOLN:
            infostring = "Maximum iterations reached with no feasible solution found";
            break;
        case MI_INFEASIBLE:
            infostring = "Problem is infeasible";
            break;
        default:
            infostring = "UNKNOWN PROBLEM IN BRANCH AND BOUND SOLVER";
    }
  } else {
    /* info dict */
    /* infostring */
    switch ( exitcode ) {
    case ECOS_OPTIMAL:
      infostring = "Optimal solution found";
      break;
    case ECOS_OPTIMAL + ECOS_INACC_OFFSET:
      infostring = "Close to optimal solution found";
      break;
    case ECOS_MAXIT:
      infostring = "Maximum number of iterations reached";
      break;
    case ECOS_PINF:
      infostring = "Primal infeasible";
      break;
    case ECOS_PINF + ECOS_INACC_OFFSET:
      infostring = "Close to primal infeasible";
      break;
    case ECOS_DINF:
      infostring = "Dual infeasible";
      break;
    case ECOS_DINF + ECOS_INACC_OFFSET:
      infostring = "Close to dual infeasible";
      break;
    case ECOS_NUMERICS:
      infostring = "Ran into numerical problems";
      break;
    case ECOS_OUTCONE:
      infostring = "PROBLEM: Multipliers leaving the cone";
      break;
    case ECOS_FATAL:
      infostring = "PROBLEM: Fatal error during initialization";
      break;
    default:
      infostring = "UNKNOWN PROBLEM IN SOLVER";
    }
  }
  numerr = 0;
  /* numerical errors */
  if ( (exitcode == ECOS_NUMERICS) || (exitcode == ECOS_OUTCONE) || (exitcode == ECOS_FATAL) ){
    numerr = 1;
  }

  /* timings */
#if PROFILING > 0
  double *tinfosVec;
  SEXP tinfos, tinfosNames;
#if PROFILING > 1
  PROTECT(tinfos = allocVector(REALSXP, 8));
  nProtected++;
  PROTECT(tinfosNames = allocVector(STRSXP, 8));
  nProtected++;
#else
  PROTECT(tinfos = allocVector(REALSXP, 3));
  nProtected++;
  PROTECT(tinfosNames = allocVector(STRSXP, 3));
  nProtected++;
#endif
  tinfosVec = REAL(tinfos);
  i = 0;
#if PROFILING > 1
  tinfosVec[i] = (double)mywork->info->tkktcreate;
  SET_STRING_ELT(tinfosNames, i, mkChar("tkktcreate"));
  i++;
  tinfosVec[i] = (double)mywork->info->tkktsolve;
  SET_STRING_ELT(tinfosNames, i, mkChar("tkktsolve"));
  i++;
  tinfosVec[i] = (double)mywork->info->tfactor;
  SET_STRING_ELT(tinfosNames, i, mkChar("tkktfactor"));
  i++;
  tinfosVec[i] = (double)mywork->info->torder;
  SET_STRING_ELT(tinfosNames, i, mkChar("torder"));
  i++;
  tinfosVec[i] = (double)mywork->info->ttranspose;
  SET_STRING_ELT(tinfosNames, i, mkChar("ttranspose"));
  i++;
#else
  tinfosVec[i] = (double)mywork->info->tsolve + (double)mywork->info->tsetup;
  SET_STRING_ELT(tinfosNames, i, mkChar("runtime"));
  i++;
  tinfosVec[i] = (double)mywork->info->tsetup;
  SET_STRING_ELT(tinfosNames, i, mkChar("tsetup"));
  i++;
  tinfosVec[i] = (double)mywork->info->tsolve;
  SET_STRING_ELT(tinfosNames, i, mkChar("tsolve"));
#endif
  setAttrib(tinfos, R_NamesSymbol, tinfosNames);
  UNPROTECT(1); /* The names have a reference now, so drop protection */
  nProtected--;
#endif

  SEXP details, detailsNames;
  double *dVec;
  PROTECT(details = allocVector(REALSXP, 11));
  nProtected++;
  PROTECT(detailsNames = allocVector(STRSXP, 11));
  nProtected++;
  dVec = REAL(details);
  dVec[0] = (double)mywork->info->pcost;
  SET_STRING_ELT(detailsNames, 0, mkChar("pcost"));
  dVec[1] = (double)mywork->info->dcost;
  SET_STRING_ELT(detailsNames, 1, mkChar("dcost"));
  dVec[2] = (double)mywork->info->pres;
  SET_STRING_ELT(detailsNames, 2, mkChar("pres"));
  dVec[3] = (double)mywork->info->dres;
  SET_STRING_ELT(detailsNames, 3, mkChar("dres"));
  dVec[4] = (double)mywork->info->pinf;
  SET_STRING_ELT(detailsNames, 4, mkChar("pinf"));
  dVec[5] = (double)mywork->info->dinf;
  SET_STRING_ELT(detailsNames, 5, mkChar("dinf"));
  dVec[6] = (double)mywork->info->pinfres;
  SET_STRING_ELT(detailsNames, 6, mkChar("pinfres"));
  dVec[7] = (double)mywork->info->dinfres;
  SET_STRING_ELT(detailsNames, 7, mkChar("dinfres"));
  dVec[8] = (double)mywork->info->gap;
  SET_STRING_ELT(detailsNames, 8, mkChar("gap"));
  dVec[9] = (double)mywork->info->relgap;
  SET_STRING_ELT(detailsNames, 9, mkChar("relgap"));
  dVec[10] = (double)mywork->stgs->feastol;
  SET_STRING_ELT(detailsNames, 10, mkChar("r0"));
  setAttrib(details, R_NamesSymbol, detailsNames);
  UNPROTECT(1); /* The names have a reference now, so drop protection */
  nProtected--;

  /* Prepare infostring for return */
  SEXP message;
  PROTECT(message = allocVector(STRSXP, 1));
  nProtected++;
  SET_STRING_ELT(message, 0, mkChar(infostring));

  /* Prepare integer vector of return codes etc. */
  SEXP intCodes, intCodesNames;
  PROTECT(intCodes = allocVector(INTSXP, 4));
  nProtected++;
  PROTECT(intCodesNames = allocVector(STRSXP, 4));
  nProtected++;
  INTEGER(intCodes)[0] = exitcode;
  SET_STRING_ELT(intCodesNames, 0, mkChar("exitFlag"));
  INTEGER(intCodes)[1] = mywork->info->iter;
  SET_STRING_ELT(intCodesNames, 1, mkChar("iter"));
  INTEGER(intCodes)[2] = mi_iterations;
  SET_STRING_ELT(intCodesNames, 2, mkChar("mi_iter"));
  INTEGER(intCodes)[3] = numerr;
  SET_STRING_ELT(intCodesNames, 3, mkChar("numerr"));
  setAttrib(intCodes, R_NamesSymbol, intCodesNames);
  UNPROTECT(1); /* The names have a reference now, so drop protection */
  nProtected--;

  /* Prepare final result to return */
  /* Number of items is number of protected items*/
  SEXP result, resultNames;
  PROTECT(result = allocVector(VECSXP, nProtected ));
  PROTECT(resultNames = allocVector(STRSXP, nProtected ));
  nProtected += 2;

  i = 0;
  SET_VECTOR_ELT(result, i, x);
  SET_STRING_ELT(resultNames, i, mkChar("x"));
  i++;
  SET_VECTOR_ELT(result, i, y);
  SET_STRING_ELT(resultNames, i, mkChar("y"));
  i++;
  SET_VECTOR_ELT(result, i, s);
  SET_STRING_ELT(resultNames, i, mkChar("s"));
  i++;
  SET_VECTOR_ELT(result, i, z);
  SET_STRING_ELT(resultNames, i, mkChar("z"));
  i++;
  SET_VECTOR_ELT(result, i, message);
  SET_STRING_ELT(resultNames, i, mkChar("infostring"));
  i++;
  SET_VECTOR_ELT(result, i, intCodes);
  SET_STRING_ELT(resultNames, i, mkChar("retcodes"));
  i++;
  SET_VECTOR_ELT(result, i, details);
  SET_STRING_ELT(resultNames, i, mkChar("summary"));
  i++;
#if PROFILING > 0
  SET_VECTOR_ELT(result, i, tinfos);
  SET_STRING_ELT(resultNames, i, mkChar("timing"));
  i++;
#endif
  setAttrib(result, R_NamesSymbol, resultNames);

  /* Clean up */
  if (num_bool > 0 || num_int > 0){
    ECOS_BB_cleanup(myecos_bb_work, 0);
  } else {
    ECOS_cleanup(mywork, 0);
  }

  free_allocated(&allocated);

  UNPROTECT(nProtected);
  return(result);

}
