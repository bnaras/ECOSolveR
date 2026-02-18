#include <R.h>
#include <Rinternals.h>

#define IGNORE_UNUSED_FUNCTIONS

#include "ecos.h"
#include "ecos_bb.h"


/** Modified from R-exts example for getting a named element */
SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  if (names == R_NilValue) return R_NilValue;
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

/* Make a copy of int vector to an idxint vector.
 * For length-0 inputs, returns a valid non-NULL pointer (allocated 1 byte)
 * so that ECOS_setup sees non-NULL and creates a valid empty spmat. */
idxint * int2idxint(SEXP intsxp) {
  int n, *source;
  idxint *result;
  n = length(intsxp);
  source = INTEGER(intsxp);
  result = (idxint *) malloc((n > 0 ? n : 1) * sizeof(idxint));
  if (result == NULL) {
    Rf_error("Out of memory in int2idxint (could not allocate %d elements)", n);
  }
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

/** Build R result list from a solved ECOS workspace.
 *  Deep-copies x/y/s/z so the workspace can be reused or freed.
 *  is_bb: nonzero if this was a branch-and-bound solve.
 *  mi_iterations: B&B iteration count (-1 for continuous). */
static SEXP build_result(pwork *mywork, idxint exitcode,
			 int is_bb, int mi_iterations) {
  const char* infostring;
  int i, nProtected = 0, numerr;
  idxint n = mywork->n, m = mywork->m, p = mywork->p;

  /* create output (all data is *deep copied*)
   * Guard: REAL() on a zero-length REALSXP may return a misaligned
   * sentinel pointer (e.g. 0x1); passing that to memcpy is UB even
   * when nbytes == 0.  Skip the copy when the dimension is zero. */
  SEXP x;
  PROTECT(x = allocVector(REALSXP, n));
  nProtected++;
  if (n > 0) memcpy(REAL(x), mywork->x, n * sizeof(double));

  SEXP y;
  PROTECT(y = allocVector(REALSXP, p));
  nProtected++;
  if (p > 0) memcpy(REAL(y), mywork->y, p * sizeof(double));

  SEXP s;
  PROTECT(s = allocVector(REALSXP, m));
  nProtected++;
  if (m > 0) memcpy(REAL(s), mywork->s, m * sizeof(double));

  SEXP z;
  PROTECT(z = allocVector(REALSXP, m));
  nProtected++;
  if (m > 0) memcpy(REAL(z), mywork->z, m * sizeof(double));

  /* Infostring based on return code */
  if (is_bb) {
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
  UNPROTECT(1);
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
  UNPROTECT(1);
  nProtected--;

  SEXP message;
  PROTECT(message = allocVector(STRSXP, 1));
  nProtected++;
  SET_STRING_ELT(message, 0, mkChar(infostring));

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
  UNPROTECT(1);
  nProtected--;

  SEXP result, resultNames;
  PROTECT(result = allocVector(VECSXP, nProtected));
  PROTECT(resultNames = allocVector(STRSXP, nProtected));
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

  UNPROTECT(nProtected);
  return result;
}

SEXP ecos_csolve(SEXP sMNP, SEXP sl, SEXP sq, SEXP se,
		 SEXP sGpr, SEXP sGjc, SEXP sGir,
		 SEXP sApr, SEXP sAjc, SEXP sAir,
		 SEXP sc, SEXP sh, SEXP sb,
		 SEXP sbool_vars, SEXP sint_vars,
		 SEXP ecosControl) {
  idxint exitcode;
  idxint n, m, p, l, ncones, *q, e, *Gjc, *Gir, *Ajc, *Air, *bool_vars, *int_vars;
  pfloat *Gpr, *Apr, *c, *h, *b;
  alloc_ptr allocated;
  int mi_iterations = -1, num_bool, num_int, is_bb;

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
  SEXP result;

  is_bb = (num_bool > 0 || num_int > 0);

  if (is_bb) {
    settings_bb opts_ecos_bb;
    setup_ECOS_BB_params(&opts_ecos_bb, ecosControl);

    myecos_bb_work = ECOS_BB_setup(n, m, p, l, ncones, q, e, Gpr, Gjc, Gir,
				   Apr, Ajc, Air, c, h, b,
				   num_bool, bool_vars, num_int, int_vars, &opts_ecos_bb);
    if ( myecos_bb_work == NULL ) {
      free_allocated(&allocated);
      Rf_error("Internal problem occurred in ECOS_BB while setting up the problem.");
      return R_NilValue; /* not reached */
    }

    mywork = myecos_bb_work->ecos_prob;
    setup_ECOS_params(mywork, ecosControl);
    exitcode = ECOS_BB_solve(myecos_bb_work);
    mi_iterations = (int) myecos_bb_work->iter;

  } else {
    mywork = ECOS_setup(n, m, p, l, ncones, q, e, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
    if ( mywork == NULL ) {
      free_allocated(&allocated);
      Rf_error("Internal problem occurred in ECOS while setting up the problem.");
      return R_NilValue; /* not reached */
    }

    setup_ECOS_params(mywork, ecosControl);
    exitcode = ECOS_solve(mywork);
  }

  result = build_result(mywork, exitcode, is_bb, mi_iterations);

  /* Clean up */
  if (is_bb) {
    ECOS_BB_cleanup(myecos_bb_work, 0);
  } else {
    ECOS_cleanup(mywork, 0);
  }

  free_allocated(&allocated);

  return result;
}

/* ------------------------------------------------------------------ */
/* Multi-step lifecycle API: setup / solve / update / cleanup         */
/* ------------------------------------------------------------------ */

/** Wrapper holding ECOS workspace and bookkeeping for the R external pointer.
 *
 *  Protected slot layout (VECSXP of length 10):
 *    [0..4]: "original" R vectors (Gpr, Apr, c, h, b) — unequilibrated,
 *            used for substituting NULL in partial updates.
 *    [5..9]: "active" duplicates (Gpr, Apr, c, h, b) — equilibrated by
 *            ECOS in-place; these are the vectors whose REAL() memory
 *            pwork borrows. Must be kept alive for the workspace lifetime.
 */
typedef struct ecos_workspace {
  pwork *work;          /* ECOS workspace (NULL after cleanup) */
  alloc_ptr allocated;  /* malloc'd idxint copies of index arrays */
  idxint n, m, p;       /* cached dimensions for update validation */
  idxint nnz_G, nnz_A;  /* cached nnz counts for update validation */
} ecos_workspace;

/** Finalizer called by R's GC when external pointer is collected. */
static void ecos_workspace_finalizer(SEXP ext_ptr) {
  ecos_workspace *ws = (ecos_workspace *)R_ExternalPtrAddr(ext_ptr);
  if (ws == NULL) return;
  if (ws->work != NULL) ECOS_cleanup(ws->work, 0);
  free_allocated(&ws->allocated);
  free(ws);
  R_ClearExternalPtr(ext_ptr);
}

/** Helper: extract ecos_workspace* from external pointer, error if NULL. */
static ecos_workspace *get_workspace(SEXP ext_ptr, const char *fname) {
  ecos_workspace *ws = (ecos_workspace *)R_ExternalPtrAddr(ext_ptr);
  if (ws == NULL || ws->work == NULL) {
    Rf_error("%s: workspace has been cleaned up", fname);
  }
  return ws;
}

/** Helper: duplicate a REAL SEXP. Returns R_NilValue if input is R_NilValue. */
static SEXP dup_real_or_nil(SEXP s) {
  if (Rf_isNull(s)) return R_NilValue;
  return Rf_duplicate(s);
}

/** Set up an ECOS workspace and return an external pointer.
 *
 *  The R-level ECOS_setup() guarantees that even when A is absent,
 *  valid empty CSC arrays (Apr=numeric(0), Ajc=integer(n+1), Air=integer(0),
 *  b=numeric(0)) are passed instead of NULL. This ensures w->A is a valid
 *  spmat (with nnz=0), which is required for ECOS_updateData's equilibration
 *  check. */
SEXP ecos_setup_R(SEXP sMNP, SEXP sl, SEXP sq, SEXP se,
		  SEXP sGpr, SEXP sGjc, SEXP sGir,
		  SEXP sApr, SEXP sAjc, SEXP sAir,
		  SEXP sc, SEXP sh, SEXP sb,
		  SEXP ecosControl) {
  idxint n, m, p, l, ncones, *q_arr, e, *Gjc, *Gir, *Ajc, *Air;
  int nProt = 0;

  ecos_workspace *ws = (ecos_workspace *)malloc(sizeof(ecos_workspace));
  if (ws == NULL) Rf_error("ecos_setup_R: out of memory");
  memset(ws, 0, sizeof(ecos_workspace));

  n = (idxint)(INTEGER(sMNP)[1]);
  m = (idxint)(INTEGER(sMNP)[0]);
  p = (idxint)(INTEGER(sMNP)[2]);

  l = (idxint)(INTEGER(sl)[0]);
  ncones = (idxint)length(sq);
  if (Rf_isNull(sq)) {
    q_arr = NULL;
  } else {
    q_arr = int2idxint(sq);
    ws->allocated.q = q_arr;
  }
  e = (idxint)(INTEGER(se)[0]);

  /* Duplicate data vectors so ECOS can equilibrate in-place
   * without modifying the caller's R objects. */
  SEXP dGpr, dApr, dc, dh, db;
  PROTECT(dGpr = dup_real_or_nil(sGpr)); nProt++;
  PROTECT(dApr = dup_real_or_nil(sApr)); nProt++;
  PROTECT(dc   = dup_real_or_nil(sc));   nProt++;
  PROTECT(dh   = dup_real_or_nil(sh));   nProt++;
  PROTECT(db   = dup_real_or_nil(sb));   nProt++;

  pfloat *Gpr = Rf_isNull(dGpr) ? NULL : (pfloat *)REAL(dGpr);
  if (Rf_isNull(sGjc)) {
    Gjc = NULL;
  } else {
    Gjc = int2idxint(sGjc);
    ws->allocated.Gjc = Gjc;
  }
  if (Rf_isNull(sGir)) {
    Gir = NULL;
  } else {
    Gir = int2idxint(sGir);
    ws->allocated.Gir = Gir;
  }

  pfloat *Apr = Rf_isNull(dApr) ? NULL : (pfloat *)REAL(dApr);
  if (Rf_isNull(sAjc)) {
    Ajc = NULL;
  } else {
    Ajc = int2idxint(sAjc);
    ws->allocated.Ajc = Ajc;
  }
  if (Rf_isNull(sAir)) {
    Air = NULL;
  } else {
    Air = int2idxint(sAir);
    ws->allocated.Air = Air;
  }

  pfloat *cptr = Rf_isNull(dc) ? NULL : (pfloat *)REAL(dc);
  pfloat *hptr = Rf_isNull(dh) ? NULL : (pfloat *)REAL(dh);
  pfloat *bptr = Rf_isNull(db) ? NULL : (pfloat *)REAL(db);

  pwork *work = ECOS_setup(n, m, p, l, ncones, q_arr, e,
			   Gpr, Gjc, Gir, Apr, Ajc, Air,
			   cptr, hptr, bptr);
  if (work == NULL) {
    free_allocated(&ws->allocated);
    free(ws);
    UNPROTECT(nProt);
    Rf_error("Internal problem occurred in ECOS while setting up the problem.");
    return R_NilValue; /* not reached */
  }

  /* Apply control settings */
  setup_ECOS_params(work, ecosControl);

  ws->work = work;
  ws->n = n;
  ws->m = m;
  ws->p = p;
  ws->nnz_G = Rf_isNull(sGpr) ? 0 : (idxint)length(sGpr);
  ws->nnz_A = Rf_isNull(sApr) ? 0 : (idxint)length(sApr);

  /* Protected slot: originals [0..4] + active duplicates [5..9] */
  SEXP prot;
  PROTECT(prot = allocVector(VECSXP, 10)); nProt++;
  SET_VECTOR_ELT(prot, 0, sGpr);   /* original Gpr */
  SET_VECTOR_ELT(prot, 1, sApr);   /* original Apr */
  SET_VECTOR_ELT(prot, 2, sc);     /* original c   */
  SET_VECTOR_ELT(prot, 3, sh);     /* original h   */
  SET_VECTOR_ELT(prot, 4, sb);     /* original b   */
  SET_VECTOR_ELT(prot, 5, dGpr);   /* active Gpr (borrowed by ECOS) */
  SET_VECTOR_ELT(prot, 6, dApr);   /* active Apr */
  SET_VECTOR_ELT(prot, 7, dc);     /* active c   */
  SET_VECTOR_ELT(prot, 8, dh);     /* active h   */
  SET_VECTOR_ELT(prot, 9, db);     /* active b   */

  SEXP ext_ptr;
  PROTECT(ext_ptr = R_MakeExternalPtr(ws, R_NilValue, prot)); nProt++;
  R_RegisterCFinalizer(ext_ptr, ecos_workspace_finalizer);

  /* Set class attribute for R dispatch */
  SEXP cls;
  PROTECT(cls = allocVector(STRSXP, 1)); nProt++;
  SET_STRING_ELT(cls, 0, mkChar("ecos_workspace"));
  setAttrib(ext_ptr, R_ClassSymbol, cls);

  UNPROTECT(nProt);
  return ext_ptr;
}

/** Solve an ECOS workspace. Optionally apply control settings first. */
SEXP ecos_solve_R(SEXP ext_ptr, SEXP ecosControl) {
  ecos_workspace *ws = get_workspace(ext_ptr, "ecos_solve_R");

  /* Optionally override settings before solve */
  if (!Rf_isNull(ecosControl)) {
    setup_ECOS_params(ws->work, ecosControl);
  }

  idxint exitcode = ECOS_solve(ws->work);
  return build_result(ws->work, exitcode, 0, -1);
}

/** Update numerical data in an ECOS workspace.
 *
 *  For any field passed as R_NilValue, the stored original (unequilibrated)
 *  data is substituted. All five data vectors are always passed to
 *  ECOS_updateData as fresh duplicates, ensuring correct equilibration
 *  behavior. */
SEXP ecos_update_R(SEXP ext_ptr, SEXP sGpr, SEXP sApr,
		   SEXP sc, SEXP sh, SEXP sb) {
  ecos_workspace *ws = get_workspace(ext_ptr, "ecos_update_R");
  SEXP prot = R_ExternalPtrProtected(ext_ptr);
  int nProt = 0;

  /* Validate lengths for user-provided (non-NULL) vectors */
  if (!Rf_isNull(sGpr) && (idxint)length(sGpr) != ws->nnz_G)
    Rf_error("ecos_update_R: Gpr length (%d) does not match original (%d)",
	     (int)length(sGpr), (int)ws->nnz_G);
  if (!Rf_isNull(sApr) && (idxint)length(sApr) != ws->nnz_A)
    Rf_error("ecos_update_R: Apr length (%d) does not match original (%d)",
	     (int)length(sApr), (int)ws->nnz_A);
  if (!Rf_isNull(sc) && (idxint)length(sc) != ws->n)
    Rf_error("ecos_update_R: c length (%d) does not match n (%d)",
	     (int)length(sc), (int)ws->n);
  if (!Rf_isNull(sh) && (idxint)length(sh) != ws->m)
    Rf_error("ecos_update_R: h length (%d) does not match m (%d)",
	     (int)length(sh), (int)ws->m);
  if (!Rf_isNull(sb) && (idxint)length(sb) != ws->p)
    Rf_error("ecos_update_R: b length (%d) does not match p (%d)",
	     (int)length(sb), (int)ws->p);

  /* For NULL args, substitute stored originals (slots 0..4) */
  SEXP src_Gpr = Rf_isNull(sGpr) ? VECTOR_ELT(prot, 0) : sGpr;
  SEXP src_Apr = Rf_isNull(sApr) ? VECTOR_ELT(prot, 1) : sApr;
  SEXP src_c   = Rf_isNull(sc)   ? VECTOR_ELT(prot, 2) : sc;
  SEXP src_h   = Rf_isNull(sh)   ? VECTOR_ELT(prot, 3) : sh;
  SEXP src_b   = Rf_isNull(sb)   ? VECTOR_ELT(prot, 4) : sb;

  /* Duplicate all five so ECOS can equilibrate fresh copies */
  SEXP dGpr, dApr, dc, dh, db;
  PROTECT(dGpr = dup_real_or_nil(src_Gpr)); nProt++;
  PROTECT(dApr = dup_real_or_nil(src_Apr)); nProt++;
  PROTECT(dc   = dup_real_or_nil(src_c));   nProt++;
  PROTECT(dh   = dup_real_or_nil(src_h));   nProt++;
  PROTECT(db   = dup_real_or_nil(src_b));   nProt++;

  pfloat *new_Gpr = Rf_isNull(dGpr) ? NULL : (pfloat *)REAL(dGpr);
  pfloat *new_Apr = Rf_isNull(dApr) ? NULL : (pfloat *)REAL(dApr);
  pfloat *new_c   = Rf_isNull(dc)   ? NULL : (pfloat *)REAL(dc);
  pfloat *new_h   = Rf_isNull(dh)   ? NULL : (pfloat *)REAL(dh);
  pfloat *new_b   = Rf_isNull(db)   ? NULL : (pfloat *)REAL(db);

  ECOS_updateData(ws->work, new_Gpr, new_Apr, new_c, new_h, new_b);

  /* Update protected slot: new originals and active duplicates */
  if (!Rf_isNull(sGpr)) SET_VECTOR_ELT(prot, 0, sGpr);
  if (!Rf_isNull(sApr)) SET_VECTOR_ELT(prot, 1, sApr);
  if (!Rf_isNull(sc))   SET_VECTOR_ELT(prot, 2, sc);
  if (!Rf_isNull(sh))   SET_VECTOR_ELT(prot, 3, sh);
  if (!Rf_isNull(sb))   SET_VECTOR_ELT(prot, 4, sb);
  SET_VECTOR_ELT(prot, 5, dGpr);
  SET_VECTOR_ELT(prot, 6, dApr);
  SET_VECTOR_ELT(prot, 7, dc);
  SET_VECTOR_ELT(prot, 8, dh);
  SET_VECTOR_ELT(prot, 9, db);

  UNPROTECT(nProt);
  return R_NilValue;
}

/** Explicitly clean up an ECOS workspace. */
SEXP ecos_cleanup_R(SEXP ext_ptr) {
  ecos_workspace *ws = (ecos_workspace *)R_ExternalPtrAddr(ext_ptr);
  if (ws == NULL) return R_NilValue; /* already cleaned up */
  if (ws->work != NULL) {
    ECOS_cleanup(ws->work, 0);
    ws->work = NULL;
  }
  free_allocated(&ws->allocated);
  free(ws);
  R_ClearExternalPtr(ext_ptr);
  return R_NilValue;
}
