context("ECOS input validation and helpers")

test_that("ecos.control returns valid defaults", {
    ctrl <- ecos.control()
    expect_type(ctrl, "list")
    expect_equal(ctrl$MAXIT, 100L)
    expect_equal(ctrl$VERBOSE, 0L)
    expect_equal(ctrl$FEASTOL, 1e-8)
    expect_equal(ctrl$ABSTOL, 1e-8)
    expect_equal(ctrl$RELTOL, 1e-8)
    expect_equal(ctrl$MI_MAX_ITERS, 1000L)
    expect_null(checkOptions(ctrl))
})

test_that("make_csc_matrix works for a dense matrix", {
    m <- matrix(c(1, 0, 0, 2, 3, 0), nrow = 2)
    csc <- make_csc_matrix(m)
    expect_type(csc, "list")
    expect_true("matbeg" %in% names(csc))
    expect_true("matind" %in% names(csc))
    expect_true("values" %in% names(csc))
    expect_equal(length(csc$values), 3L)
})

test_that("mismatched dimensions error correctly", {
    if (!requireNamespace("Matrix", quietly = TRUE)) skip("Matrix not available")
    ## G has 3 rows, but l=1, q=NULL, e=0 sums to 1, not 3
    G <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L), x = c(-1, -1), dims = c(3L, 2L))
    cc <- c(1.0, 1.0)
    h <- c(0.0, 0.0, 0.0)
    dims <- list(l = 1L, q = NULL, e = 0L)
    expect_error(ECOS_csolve(c = cc, G = G, h = h, dims = dims),
                 "does not match")
})

test_that("plain matrix A works when G is dgCMatrix", {
    if (!requireNamespace("Matrix", quietly = TRUE)) skip("Matrix not available")
    ## Simple LP: minimize -x1 - x2 s.t. x1 + x2 = 1, x1 >= 0, x2 >= 0
    G <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L), x = c(-1, -1), dims = c(2L, 2L))
    cc <- c(-1.0, -1.0)
    h <- c(0.0, 0.0)
    dims <- list(l = 2L, q = NULL, e = 0L)
    A <- matrix(c(1, 1), nrow = 1)
    b <- 1.0
    retval <- ECOS_csolve(c = cc, G = G, h = h, dims = dims, A = A, b = b)
    expect_false(is.null(retval))
    expect_equal(retval$retcodes[["exitFlag"]], 0L)
    expect_equal(sum(retval$x), 1.0, tolerance = 1e-6)
})
