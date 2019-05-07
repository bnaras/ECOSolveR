context("ECOS NAN Tests")
test_that("Test that convolution example works on all platforms including 32-bit", {
    n <- 5L
    c <- rep(-2.271944, n)
    dims <- list(f = 0L, l = 0L, q = NULL, e = 0L)
    G <- Matrix::Matrix(numeric(0), nrow = 0, ncol = n, sparse = TRUE)
    A <- Matrix::Matrix(numeric(0), nrow = 0, ncol = n, sparse = TRUE)
    retval <- ECOS_csolve(c = c,
                          G = G, h = numeric(0),
                          dims = dims,
                          A = A, b = numeric(0))
    ## Answer: Dual infeasible
    expect_equal(c(exitFlag = 2L), retval$retcodes["exitFlag"])
})
