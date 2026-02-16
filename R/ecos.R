#' Convert a plain matrix or simple triplet form matrix to a [Matrix::dgCMatrix-class] (implicit) form
#' @param x a matrix or a simple triplet form matrix
#' @return a list of row pointer, column pointer, and values corresponding to a [Matrix::dgCMatrix-class] object
#' @keywords internal
make_csc_matrix <- function(x) UseMethod("make_csc_matrix")

#' @method make_csc_matrix matrix
make_csc_matrix.matrix <- function(x) {
    if( !is.matrix(x) )
        cli_abort("{.arg x} must be a matrix.")

    ind <- which(x != 0, arr.ind = TRUE)
    list(matbeg = c(0L, cumsum(tabulate(ind[, 2L], ncol(x)))),
         matind = ind[, 1L] - 1L,
         values = x[ind])
}

#' @method make_csc_matrix simple_triplet_matrix
make_csc_matrix.simple_triplet_matrix <- function(x) {
    if(!inherits(x, "simple_triplet_matrix"))
        cli_abort("{.arg x} must be of class {.cls simple_triplet_matrix}.")

    ## The matrix method assumes that indices for non-zero entries are
    ## in row-major order, but the simple_triplet_matrix() constructor
    ## currently does not canonicalize accordingly ...
    ind <- order(x$j, x$i)
    list(matbeg = c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
         matind = x$i[ind] - 1L,
         values = x$v[ind])
}

## Validate and prepare problem data for ECOS C calls.
## Returns a list with C-ready components:
##   MNP (integer[3]), l, q, e, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b
.validate_and_prepare <- function(c, G, h, dims, A, b) {
    dims <- lapply(dims, as.integer)

    nullG <- (is.null(G) || prod(dim(G)) == 0L)
    nontrivialH <- isNontrivialNumericVector(h)

    if ((nullG && nontrivialH) ||
        (!nullG && !nontrivialH)) {
        cli_abort("{.arg G} and {.arg h} must be supplied together.")
    }

    nullA <- (is.null(A) || prod(dim(A)) == 0L)
    nontrivialB <- isNontrivialNumericVector(b)
    if ((nullA && nontrivialB) ||
        (!nullA && !nontrivialB)) {
        cli_abort("{.arg A} and {.arg b} must be supplied together.")
    }

    nC <- length(c)

    if (nullG) {
        Gpr <- h <- numeric(0)
        Gir <- integer(0)
        Gjc <- integer(nC + 1L)
        mG <- 0L
        nG <- nC
    } else {
        if (inherits(G, c("matrix", "simple_triplet_matrix"))) {
            csc <- make_csc_matrix(G)
            Gpr <- csc[["values"]]
            Gir <- csc[["matind"]]
            Gjc <- csc[["matbeg"]]
        } else if (inherits(G, "dgCMatrix")) {
            Gpr <- G@x
            Gir <- G@i
            Gjc <- G@p
        } else {
            cli_abort("{.arg G} must be a {.cls dgCMatrix}, {.cls matrix}, or {.cls simple_triplet_matrix}.")
        }
        mG <- nrow(G)
        nG <- ncol(G)
        if (nG != nC) {
            cli_abort("Number of columns of {.arg G} must match length of {.arg c}.")
        }
    }

    if (nullA) {
        Apr <- Air <- Ajc <- b <- NULL
        mA <- nA <- 0L
    } else {
        if (inherits(A, "sparseMatrix")) {
            if (!inherits(A, "dgCMatrix")) A  <- as(as(A, "CsparseMatrix"), "dgCMatrix")
            Apr <- A@x
            Air <- A@i
            Ajc <- A@p
        } else if (inherits(A, c("matrix", "simple_triplet_matrix"))) {
            csc <- make_csc_matrix(A)
            Apr <- csc[["values"]]
            Air <- csc[["matind"]]
            Ajc <- csc[["matbeg"]]
        } else {
            cli_abort("{.arg A} must be a {.cls dgCMatrix}, {.cls matrix}, or {.cls simple_triplet_matrix}.")
        }
        mA <- nrow(A)
        nA <- ncol(A)
        if (mA != length(b)) {
            cli_abort("{.arg b} has incompatible dimension with {.arg A}.")
        }
        if (nA != nC) {
            cli_abort("Number of columns of {.arg A} must match length of {.arg c}.")
        }
    }

    ## Need to check dims as well
    if (is.null(dims)) {
        cli_abort("{.arg dims} must be a non-null list.")
    }
    ## dimensions of the positive orthant cone
    l <- dims$l
    if (is.null(l)) {
        l <- 0L
    } else {
        if (!isNonnegativeInt(l))
            cli_abort("{.code dims$l} must be a non-negative integer.")
    }
    ## dimensions of the second order cones
    q <- dims$q
    if (!is.null(q)) {
        if (typeof(q) != "integer" || !all(q > 0L))
            cli_abort("{.code dims$q} must be an integer vector of positive values.")
    }
    ## number of exponential cones
    e <- dims$e
    if (is.null(e)) {
        e <- 0L
    } else {
        if (!isNonnegativeInt(e))
            cli_abort("{.code dims$e} must be a non-negative integer.")
    }
    ## check that sum(q) + l + 3 * e = m
    if ( (sum(q) + l + 3L * e) != mG ) {
        cli_abort("Number of rows of {.arg G} does not match {.code dims$l + sum(dims$q) + 3 * dims$e}.")
    }

    list(MNP = c(mG, nC, mA),
         l = l, q = q, e = e,
         Gpr = Gpr, Gjc = Gjc, Gir = Gir,
         Apr = Apr, Ajc = Ajc, Air = Air,
         c = c, h = h, b = b)
}

#' Solve a conic optimization problem
#'
#' The function \code{ECOS_csolve} is a wrapper around the ecos
#' \code{csolve} C function. Conic constraints are specified using the
#' \eqn{G} and \eqn{h} parameters and can be \code{NULL} and zero
#' length vector respectively indicating an absence of conic
#' constraints.  Similarly, equality constraints are specified via
#' \eqn{A} and \eqn{b} parameters with \code{NULL} and empty vector
#' values representing a lack of such constraints. At most one of the
#' pair \eqn{(G , h)} or \eqn{(A, b)} is allowed to be absent.
#'
#' @param c the coefficients of the objective function; the length of
#'     this determines the number of variables \eqn{n} in the problem.
#' @param G the inequality constraint matrix in one of three forms: a
#'     plain matrix, simple triplet matrix, or compressed column
#'     format, e.g. \link[Matrix]{dgCMatrix-class}. Can also be
#'     \code{NULL}
#' @param h the right hand side of the inequality constraint. Can be
#'     empty numeric vector.
#' @param dims is a list of three named elements: \code{dims['l']} an
#'     integer specifying the dimension of positive orthant cone,
#'     \code{dims['q']} an integer vector specifying dimensions of
#'     second-order cones, \code{dims['e']} an integer specifying the
#'     number of exponential cones
#' @param A the optional equality constraint matrix in one of three
#'     forms: a plain matrix, simple triplet matrix, or compressed
#'     column format, e.g. \link[Matrix]{dgCMatrix-class}. Can be
#'     \code{NULL}
#' @param b the right hand side of the equality constraint, must be
#'     specified if \eqn{A} is. Can be empty numeric vector.
#' @param bool_vars the indices of the variables, 1 through \eqn{n},
#'     that are boolean; that is, they are either present or absent in
#'     the solution
#' @param int_vars the indices of the variables, 1 through \eqn{n},
#'     that are integers
#' @param control is a named list that controls various optimization
#'     parameters; see \link[ECOSolveR]{ecos.control}.
#'
#' @return a list of 8 named items
#'  \describe{
#'   \item{x}{primal variables}
#'   \item{y}{dual variables for equality constraints}
#'   \item{s}{slacks for \eqn{Gx + s <= h}, \eqn{s \in K}}
#'   \item{z}{dual variables for inequality constraints \eqn{s \in K}}
#'   \item{infostring}{gives information about the status of solution}
#'   \item{retcodes}{a named integer vector containing four elements
#'     \describe{
#'       \item{exitflag}{0=\code{ECOS_OPTIMAL}, 1=\code{ECOS_PINF},
#'          2=\code{ECOS_DINF}, 10=\code{ECOS_INACC_OFFSET}, -1=\code{ECOS_MAXIT},
#'          -2=\code{ECOS_NUMERICS}, -3=\code{ECOS_OUTCONE}, -4=\code{ECOS_SIGINT},
#'          -7=\code{ECOS_FATAL}. See \link[ECOSolveR]{ECOS_exitcodes}}.
#'       \item{iter}{the number of iterations used}
#'       \item{mi_iter}{the number of iterations for mixed integer problems}
#'       \item{numerr}{a non-zero number if a numeric error occurred}
#'     }
#'   }
#'   \item{summary}{a named numeric vector containing
#'     \describe{
#'       \item{pcost}{value of primal objective}
#'       \item{dcost}{value of dual objective}
#'       \item{pres}{primal residual on inequalities and equalities}
#'       \item{dres}{dual residual}
#'       \item{pinf}{primal infeasibility measure}
#'       \item{dinf}{dual infeasibility measure}
#'       \item{pinfres}{primal infeasibility residual}
#'       \item{dinfres}{dual infeasibility residual}
#'       \item{gap}{duality gap}
#'       \item{relgap}{relative duality gap}
#'       \item{r0}{Unknown at the moment to this R package maintainer.}
#'     }
#'   }
#'   \item{timing}{a named numeric vector of timing information consisting of
#'     \describe{
#'       \item{runtime}{the total runtime in ecos}
#'       \item{tsetup}{the time for setup of the problem}
#'       \item{tsolve}{the time to solve the problem}
#'     }
#'   }
#' }
#'
#' @section Details:
#'
#' A call to this function will solve the problem:
#' minimize \eqn{c^Tx}, subject to \eqn{Ax = b}, and \eqn{h - G*x \in K}.
#'
#' Variables can be constrained to be boolean (1 or 0) or integers. This is indicated
#' by specifying parameters \code{bool_vars} and/or \code{int_vars} respectively. If so
#' indicated, the solutions will be found using a branch and bound algorithm.
#'
#' @examples
#'
#' ## githubIssue98
#' cat("Basic matrix interface\n")
#' Gmat <- matrix(c(0.416757847405471, 2.13619609566845, 1.79343558519486, 0, 0,
#'                  0, 0, -1, 0, 0, 0, 0.056266827226329, -1.64027080840499, 0.841747365656204,
#'                  0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0.416757847405471, 2.13619609566845,
#'                  1.79343558519486, 0, 0, 0, -1, 0, 0, 0, 0, 0.056266827226329, -1.64027080840499,
#'                  0.841747365656204, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0), ncol = 5L)
#' c <- as.numeric(c(0, 0, 0, 0, 1))
#' h <- as.numeric(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#' dims <- list(l = 6L, q = 5L, e = 0L)
#' ECOS_csolve(c = c, G = Gmat, h = h,
#'            dims = dims,
#'            A = NULL, b = numeric(0))
#'
#' cat("Simple Triplet Matrix interface, if you have package slam\n")
#' if (requireNamespace("slam")) {
#'
#'   ECOS_csolve(c = c, G = slam::as.simple_triplet_matrix(Gmat), h = h,
#'               dims = dims,
#'               A = NULL, b = numeric(0))
#' }
#'
#' if (requireNamespace("Matrix")) {
#'    ECOS_csolve(c = c, G = Matrix::Matrix(Gmat), h = h,
#'                dims = dims,
#'                A = NULL, b = numeric(0))
#' }
#'
#' ## Larger problems using saved data can be found in the test suite.
#' ## Here is one
#' if (requireNamespace("Matrix")) {
#'   MPC01 <- readRDS(system.file("testdata", "MPC01_1.RDS", package = "ECOSolveR"))
#'   G <- Matrix::sparseMatrix(x = MPC01$Gpr, i = MPC01$Gir, p = MPC01$Gjc,
#'                             dims = c(MPC01$m, MPC01$n), index1 = FALSE)
#'   h <- MPC01$h
#'   dims <- lapply(list(l = MPC01$l, q=MPC01$q, e=MPC01$e), as.integer)
#'   retval <- ECOS_csolve(c = MPC01$c, G=G, h = h, dims = dims, A = NULL, b = NULL,
#'                         control = ecos.control(verbose=1L))
#'   retval$retcodes
#'   retval$infostring
#'   retval$summary
#' }
#' @importFrom cli cli_abort
#' @importFrom methods as
#'
#' @export
ECOS_csolve <- function(c = numeric(0), G = NULL, h=numeric(0),
                         dims=list(l = integer(0), q = NULL, e = integer(0)),
                         A = NULL, b = numeric(0),
                         bool_vars = integer(0), int_vars = integer(0),
                         control = ecos.control() ) {

    checkOptions(control)
    prep <- .validate_and_prepare(c, G, h, dims, A, b)
    nC <- length(c)

    bool_vars <- as.integer(bool_vars)
    if (( length(bool_vars) > 0L) && any(bool_vars < 1L | bool_vars > nC) ) {
        cli_abort("{.arg bool_vars} must be integers between 1 and {nC}.")
    } else {
        bool_vars <- sort.int(bool_vars - 1L)
    }

    int_vars <- as.integer(int_vars)
    if (( length(int_vars) > 0L) && any(int_vars < 1L | int_vars > nC) ) {
        cli_abort("{.arg int_vars} must be integers between 1 and {nC}.")
    } else {
        int_vars <- sort.int(int_vars - 1L)
    }

    result <- .Call('ecos_csolve',
                    prep$MNP,
                    prep$l, prep$q, prep$e,
                    prep$Gpr, prep$Gjc, prep$Gir,
                    prep$Apr, prep$Ajc, prep$Air,
                    prep$c, prep$h, prep$b,
                    bool_vars, int_vars,
                    control,
                    PACKAGE = 'ECOSolveR')
    result
}

## ------------------------------------------------------------------ ##
## Multi-step lifecycle: ECOS_setup / ECOS_solve / ECOS_update /      ##
## ECOS_cleanup                                                       ##
## ------------------------------------------------------------------ ##

#' Set up an ECOS workspace for multi-step solving
#'
#' Creates an ECOS workspace that can be solved, updated with new
#' numerical data, and solved again without repeating the expensive
#' symbolic analysis phase. This is useful for parametric optimization,
#' model predictive control, and other settings where the problem
#' structure stays the same but data changes.
#'
#' @inheritParams ECOS_csolve
#' @return an external pointer of class \code{"ecos_workspace"}.
#'   Must eventually be freed via \code{\link{ECOS_cleanup}} or R
#'   garbage collection.
#'
#' @seealso \code{\link{ECOS_solve}}, \code{\link{ECOS_update}},
#'   \code{\link{ECOS_cleanup}}, \code{\link{ECOS_csolve}}
#' @export
ECOS_setup <- function(c, G, h,
                       dims = list(l = integer(0), q = NULL, e = integer(0)),
                       A = NULL, b = numeric(0),
                       control = ecos.control()) {
    checkOptions(control)
    prep <- .validate_and_prepare(c, G, h, dims, A, b)
    nC <- length(c)
    ## Ensure valid empty CSC arrays when A is absent.
    ## ECOS_updateData's equilibration check dereferences w->A
    ## unconditionally, so we need w->A to be a valid spmat (not NULL).
    if (is.null(prep$Apr)) {
        prep$Apr <- numeric(0)
        prep$Ajc <- integer(nC + 1L)
        prep$Air <- integer(0)
        prep$b <- numeric(0)
    }
    .Call('ecos_setup_R',
          prep$MNP,
          prep$l, prep$q, prep$e,
          prep$Gpr, prep$Gjc, prep$Gir,
          prep$Apr, prep$Ajc, prep$Air,
          prep$c, prep$h, prep$b,
          control,
          PACKAGE = 'ECOSolveR')
}

#' Solve an ECOS workspace
#'
#' Calls \code{ECOS_solve} on an existing workspace created by
#' \code{\link{ECOS_setup}}. The workspace can be solved multiple
#' times, optionally with \code{\link{ECOS_update}} calls in between
#' to change numerical data.
#'
#' @param workspace an external pointer of class \code{"ecos_workspace"}
#'   as returned by \code{\link{ECOS_setup}}.
#' @param control optional named list of solver parameters (see
#'   \code{\link{ecos.control}}). If \code{NULL} (the default), the
#'   settings from \code{\link{ECOS_setup}} (or the most recent
#'   \code{ECOS_solve} call with non-NULL control) are reused.
#' @return the same result list as \code{\link{ECOS_csolve}}.
#'
#' @seealso \code{\link{ECOS_setup}}, \code{\link{ECOS_update}},
#'   \code{\link{ECOS_cleanup}}
#' @export
ECOS_solve <- function(workspace, control = NULL) {
    if (!inherits(workspace, "ecos_workspace"))
        cli_abort("{.arg workspace} must be an {.cls ecos_workspace} object.")
    if (!is.null(control)) checkOptions(control)
    .Call('ecos_solve_R', workspace, control, PACKAGE = 'ECOSolveR')
}

#' Update numerical data in an ECOS workspace
#'
#' Replaces numerical values in an existing ECOS workspace without
#' repeating symbolic analysis. The sparsity structure must remain
#' the same; only the non-zero values of G and A, and the dense
#' vectors c, h, b can be updated. Pass \code{NULL} for any argument
#' to leave it unchanged.
#'
#' @param workspace an external pointer of class \code{"ecos_workspace"}
#'   as returned by \code{\link{ECOS_setup}}.
#' @param Gpr numeric vector of new non-zero values for G (same length
#'   as original), or \code{NULL} to keep current values.
#' @param Apr numeric vector of new non-zero values for A (same length
#'   as original), or \code{NULL} to keep current values.
#' @param c numeric vector of new objective coefficients, or \code{NULL}.
#' @param h numeric vector of new inequality RHS, or \code{NULL}.
#' @param b numeric vector of new equality RHS, or \code{NULL}.
#' @return the \code{workspace} object, invisibly.
#'
#' @seealso \code{\link{ECOS_setup}}, \code{\link{ECOS_solve}},
#'   \code{\link{ECOS_cleanup}}
#' @export
ECOS_update <- function(workspace, Gpr = NULL, Apr = NULL,
                        c = NULL, h = NULL, b = NULL) {
    if (!inherits(workspace, "ecos_workspace"))
        cli_abort("{.arg workspace} must be an {.cls ecos_workspace} object.")
    .Call('ecos_update_R', workspace, Gpr, Apr, c, h, b,
          PACKAGE = 'ECOSolveR')
    invisible(workspace)
}

#' Clean up an ECOS workspace
#'
#' Frees the C-level ECOS workspace. After this call the workspace
#' external pointer is invalidated and cannot be used for further
#' solves or updates. It is safe to call this more than once.
#'
#' @param workspace an external pointer of class \code{"ecos_workspace"}
#'   as returned by \code{\link{ECOS_setup}}.
#' @return \code{NULL}, invisibly.
#'
#' @seealso \code{\link{ECOS_setup}}, \code{\link{ECOS_solve}},
#'   \code{\link{ECOS_update}}
#' @export
ECOS_cleanup <- function(workspace) {
    .Call('ecos_cleanup_R', workspace, PACKAGE = 'ECOSolveR')
    invisible(NULL)
}

#' Print method for ECOS workspace objects
#'
#' @param x an \code{"ecos_workspace"} external pointer.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.ecos_workspace <- function(x, ...) {
    cat("<ecos_workspace>\n")
    invisible(x)
}
