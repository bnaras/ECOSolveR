context("ECOS Lifecycle API")

## Helper: load a test problem and construct sparse G matrix
load_problem <- function(name) {
    d <- readRDS(system.file("testdata", name, package = "ECOSolveR"))
    G <- Matrix::sparseMatrix(x = d$Gpr, i = d$Gir, p = d$Gjc,
                              dims = c(d$m, d$n), index1 = FALSE)
    dims <- lapply(list(l = d$l, q = d$q, e = d$e), as.integer)
    if (!is.null(d$Apr) && length(d$Apr) > 0) {
        A <- Matrix::sparseMatrix(x = d$Apr, i = d$Air, p = d$Ajc,
                                  dims = c(d$p, d$n), index1 = FALSE)
        b <- d$b
    } else {
        A <- NULL
        b <- numeric(0)
    }
    list(c = d$c, G = G, h = d$h, dims = dims, A = A, b = b, raw = d)
}

## -----------------------------------------------------------------
## 1. Basic lifecycle: setup -> solve -> cleanup
## -----------------------------------------------------------------
test_that("basic lifecycle works", {
    prob <- load_problem("update_data_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims, A = prob$A, b = prob$b)
    on.exit(ECOS_cleanup(ws), add = TRUE)
    expect_s3_class(ws, "ecos_workspace")

    res <- ECOS_solve(ws)
    expect_equal(res$retcodes[["exitFlag"]], 0L)
    expect_equal(res$infostring, "Optimal solution found")
})

## -----------------------------------------------------------------
## 2. Results match ECOS_csolve
## -----------------------------------------------------------------
test_that("lifecycle results match ECOS_csolve", {
    prob <- load_problem("update_data_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims, A = prob$A, b = prob$b)
    on.exit(ECOS_cleanup(ws), add = TRUE)

    res_lc <- ECOS_solve(ws)
    res_ref <- ECOS_csolve(c = prob$c, G = prob$G, h = prob$h,
                           dims = prob$dims, A = prob$A, b = prob$b)

    expect_equal(res_lc$x, res_ref$x, tolerance = 1e-8)
    expect_equal(res_lc$y, res_ref$y, tolerance = 1e-8)
    expect_equal(res_lc$s, res_ref$s, tolerance = 1e-8)
    expect_equal(res_lc$z, res_ref$z, tolerance = 1e-8)
    expect_equal(res_lc$retcodes, res_ref$retcodes)
})

## -----------------------------------------------------------------
## 3. Update + re-solve produces different solution
## -----------------------------------------------------------------
test_that("update and re-solve works", {
    prob <- load_problem("update_data_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims, A = prob$A, b = prob$b)
    on.exit(ECOS_cleanup(ws), add = TRUE)

    res1 <- ECOS_solve(ws)
    expect_equal(res1$retcodes[["exitFlag"]], 0L)

    ## Perturb h and update
    new_h <- prob$h + 0.5
    ECOS_update(ws, h = new_h)
    res2 <- ECOS_solve(ws)
    expect_equal(res2$retcodes[["exitFlag"]], 0L)

    ## Solutions should differ
    expect_false(isTRUE(all.equal(res1$x, res2$x)))

    ## Compare with a fresh ECOS_csolve using perturbed h
    res_ref <- ECOS_csolve(c = prob$c, G = prob$G, h = new_h,
                           dims = prob$dims, A = prob$A, b = prob$b)
    expect_equal(res2$x, res_ref$x, tolerance = 1e-6)
})

## -----------------------------------------------------------------
## 4. Wrong-length update errors
## -----------------------------------------------------------------
test_that("wrong-length update vectors error", {
    prob <- load_problem("update_data_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims, A = prob$A, b = prob$b)
    on.exit(ECOS_cleanup(ws), add = TRUE)

    expect_error(ECOS_update(ws, h = numeric(1)), "h length")
    expect_error(ECOS_update(ws, c = numeric(1)), "c length")
    expect_error(ECOS_update(ws, b = numeric(1)), "b length")
    expect_error(ECOS_update(ws, Gpr = numeric(1)), "Gpr length")
    expect_error(ECOS_update(ws, Apr = numeric(1)), "Apr length")
})

## -----------------------------------------------------------------
## 5. Use-after-cleanup errors
## -----------------------------------------------------------------
test_that("use after cleanup errors", {
    prob <- load_problem("update_data_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims, A = prob$A, b = prob$b)
    ECOS_cleanup(ws)

    expect_error(ECOS_solve(ws), "cleaned up")
    expect_error(ECOS_update(ws, h = prob$h), "cleaned up")
    ## Double cleanup should be safe (no error)
    expect_silent(ECOS_cleanup(ws))
})

## -----------------------------------------------------------------
## 6. GC finalizer doesn't crash
## -----------------------------------------------------------------
test_that("GC finalizer does not crash", {
    prob <- load_problem("update_data_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims, A = prob$A, b = prob$b)
    rm(ws)
    expect_silent(gc())
})

## -----------------------------------------------------------------
## 7. Per-solve control override
## -----------------------------------------------------------------
test_that("per-solve control override works", {
    prob <- load_problem("update_data_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims, A = prob$A, b = prob$b)
    on.exit(ECOS_cleanup(ws), add = TRUE)

    ## Solve with maxit=1 -> should hit MAXIT
    res_maxit <- ECOS_solve(ws, control = ecos.control(maxit = 1L))
    expect_equal(res_maxit$retcodes[["exitFlag"]], -1L)

    ## Re-solve with default settings -> should find optimal
    res_opt <- ECOS_solve(ws, control = ecos.control())
    expect_equal(res_opt$retcodes[["exitFlag"]], 0L)
})

## -----------------------------------------------------------------
## 8. Partial updates (update only c, only h)
## -----------------------------------------------------------------
test_that("partial updates work correctly", {
    prob <- load_problem("update_data_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims, A = prob$A, b = prob$b)
    on.exit(ECOS_cleanup(ws), add = TRUE)

    res_orig <- ECOS_solve(ws)

    ## Update only c
    new_c <- prob$c * 2
    ECOS_update(ws, c = new_c)
    res_c <- ECOS_solve(ws)
    ## Optimal value should roughly double since c is scaled
    expect_false(isTRUE(all.equal(res_orig$summary[["pcost"]],
                                  res_c$summary[["pcost"]])))

    ## Update only h back to original, keeping scaled c
    ECOS_update(ws, h = prob$h)
    res_h <- ECOS_solve(ws)
    expect_equal(res_h$retcodes[["exitFlag"]], 0L)
})

## -----------------------------------------------------------------
## 9. Problem without A (no equality constraints)
## -----------------------------------------------------------------
test_that("lifecycle works without equality constraints", {
    prob <- load_problem("MPC01_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims)
    on.exit(ECOS_cleanup(ws), add = TRUE)

    res <- ECOS_solve(ws)
    expect_equal(res$retcodes[["exitFlag"]], 0L)

    ref <- ECOS_csolve(c = prob$c, G = prob$G, h = prob$h,
                       dims = prob$dims)
    expect_equal(res$x, ref$x, tolerance = 1e-6)

    ## Update h and re-solve
    new_h <- prob$h + 0.01
    ECOS_update(ws, h = new_h)
    res2 <- ECOS_solve(ws)
    expect_false(isTRUE(all.equal(res$x, res2$x)))
})

## -----------------------------------------------------------------
## 10. Print method
## -----------------------------------------------------------------
test_that("print method works", {
    prob <- load_problem("MPC01_1.RDS")
    ws <- ECOS_setup(c = prob$c, G = prob$G, h = prob$h,
                     dims = prob$dims)
    on.exit(ECOS_cleanup(ws), add = TRUE)
    expect_output(print(ws), "ecos_workspace")
})
