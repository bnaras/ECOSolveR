
make_csc_matrix <- function(x) UseMethod("make_csc_matrix")

make_csc_matrix.matrix <- function(x) {
    if( !is.matrix(x) )
        stop("Argument 'x' must be a matrix.")

    ind <- which(x != 0, arr.ind = TRUE)
    list(matbeg = c(0L, cumsum(tabulate(ind[, 2L], ncol(x)))),
         matind = ind[, 1] - 1L,
         values = x[ind])
}

make_csc_matrix.simple_triplet_matrix <- function(x) {
    if(!inherits(x, "simple_triplet_matrix"))
        stop("Argument 'x' must be of class 'simple_triplet_matrix'.")

    ## The matrix method assumes that indices for non-zero entries are
    ## in row-major order, but the simple_triplet_matrix() constructor
    ## currently does not canonicalize accordingly ...
    ind <- order(x$j, x$i)
    list(matbeg = c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
         matind = x$i[ind] - 1L,
         values = x$v[ind])
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
#' @param h the right hand size of the inequality constraint. Can be
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
#'
#' @export
ECOS_csolve <- function(c = numeric(0), G = NULL, h=numeric(0),
                         dims=list(l = integer(0), q = NULL, e = integer(0)),
                         A = NULL, b = numeric(0),
                         bool_vars = integer(0), int_vars = integer(0),
                         control = ecos.control() ) {

    if (!is.null(optionCheck <- checkOptions(control))) {
        stop(optionCheck)
    }
    dims <- lapply(dims, as.integer)

    nullG <- (is.null(G) || prod(dim(G)) == 0L)
    nontrivialH <- isNontrivialNumericVector(h)

    if ((nullG && nontrivialH) ||
        (!nullG && !nontrivialH)) {
        stop("G and h must be supplied together")
    }

    nullA <- (is.null(A) || prod(dim(A)) == 0L)
    nontrivialB <- isNontrivialNumericVector(b)
    if ((nullA && nontrivialB) ||
        (!nullA && !nontrivialB)) {
        stop("A and b must be supplied together")
    }

    nC <- length(c)

    if (nullG) {
        Gpr <- h <- numeric(0)
        Gir <- integer(0)
        Gjc <- integer(nC + 1L)
        mG <- 0L
        nG <- nC
    } else {
        if ( inherits(G, "CsparseMatrix") ) {
            Gpr <- G@x
            Gir <- G@i
            Gjc <- G@p
        } else if (inherits(G, c("matrix", "simple_triplet_matrix"))) {
            csc <- make_csc_matrix(G)
            Gpr <- csc[["values"]]
            Gir <- csc[["matind"]]
            Gjc <- csc[["matbeg"]]
        } else {
            stop("G is required to be of class dgCMatrix or matrix or simple_triplet_matrix")
        }
        mG <- nrow(G)
        nG <- ncol(G)
        if (nG != nC) {
            stop("Column length of G and length of c should match")
        }
    }

    if (nullA) {
        Apr <- Air <- Ajc <- b <- NULL
        mA <- nA <- 0L
    } else {
        if ( inherits(G, "CsparseMatrix") ) {
            Apr <- A@x
            Air <- A@i
            Ajc <- A@p
        } else if (inherits(G, c("matrix", "simple_triplet_matrix"))) {
            csc <- make_csc_matrix(A)
            Apr <- csc[["values"]]
            Air <- csc[["matind"]]
            Ajc <- csc[["matbeg"]]
        } else {
            stop("A is required to be of class dgCMatrix")
        }
        mA <- nrow(A)
        nA <- ncol(A)
        if (mA != length(b)) {
            stop("b has incompatible dimension with A")
        }
        if (nA != nC) {
            stop("Column length of A and length of c should match")
        }
    }


    ## Need to check dims as well
    if (is.null(dims)) {
        stop("dims must be a non-null list")
    }
    ## dimensions of the positive orthant cone
    l <- dims$l
    if (is.null(l)) {
        l <- 0L
    } else {
        if (!isNonnegativeInt(l))
            stop("dims['l'] should be a non-negative int")
    }
    ## dimensions of the second order cones
    q <- dims$q
    if (!is.null(q)) {
        if (typeof(q) != "integer" || !all(q > 0L))
            stop("dims['q'] should be an integer vector of positive integers")
    }
    ## number of exponential cones
    e <- dims$e
    if (is.null(e)) {
        e <- 0L
    } else {
        if (!isNonnegativeInt(e))
            stop("dims['e'] should be a non-negative int")
    }
    ## I am not performing this check for now...
    ## check that sum(q) + l + 3 * e = m
    ## if ( (sum(q) + l + 3 * e) != m ) {
    ##     stop("Number of rows of G does not match dims$l + sum(dims$q) + dims$e");
    ## }

    bool_vars <- as.integer(bool_vars)
    if (( length(bool_vars) > 0L) && any(bool_vars < 1L | bool_vars > nC) ) {
        stop(sprintf("bool_vars must be integers between 1 and %d", nC))
    } else {
        bool_vars <- sort.int(bool_vars - 1L)
    }

    int_vars <- as.integer(int_vars)
    if (( length(int_vars) > 0L) && any(int_vars < 1L | int_vars > nC) ) {
        stop(sprintf("int_vars must be integers between 1 and %d", nC))
    } else {
        int_vars <- sort.int(int_vars - 1L)
    }

    result <- .Call('ecos_csolve',
                    c(mG, nC, mA),
                    l, q, e,
                    Gpr, Gjc, Gir,
                    Apr, Ajc, Air,
                    c, h, b,
                    bool_vars, int_vars,
                    control,
                    PACKAGE = 'ECOSolveR')
    result
}
