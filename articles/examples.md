# ECOSolveR Examples

## $`L_1`$ minimization (Linear Programming)

We solve the following problem that arises for example in sparse signal
reconstruction problems such as compressed sensing: \$\$ \mbox{minimize
} \|\|x\|\|\_1 \mbox{ (\$L_1\$) }\\ \mbox{subject to } Ax = b \$\$

with $`x\in R^n`$, $`A \in R^{m \times n}`$ and $`m\leq n.`$ Reformulate
the problem expressing the $`L_1`$ norm of $`x`$ as follows \$\$ x \leq
u\\ -x \leq u\\ \$\$

where $`u\in R^n`$ and we minimize the sum of $`u`$. The reformulated
problem using the stacked variables

``` math
z = \begin{pmatrix}x\\u\end{pmatrix}
```

is now \$\$ \mbox{minimize } c^{\top}z\\ \mbox{subject to } \tilde{A}x =
b \mbox{ (LP) }\\ Gx \leq h \$\$ where the inequality is with respective
to the positive orthant.

Here is the R code that generates a random instance of this problem and
solves it.

``` r

library(ECOSolveR)
library(Matrix)
set.seed(182391)
n <- 1000L
m <- 10L
density <- 0.01
c <- c(rep(0.0, n), rep(1.0, n))
```

First, a function to generate random sparse matrices with normal
entries.

``` r

sprandn <- function(nrow, ncol, density) {
    items <- ceiling(nrow * ncol * density)
    matrix(c(rnorm(items),
             rep(0, nrow * ncol - items)),
           nrow = nrow)
}
```

``` r

A <- sprandn(m, n, density)
Atilde <- Matrix(cbind(A, matrix(rep(0.0, m * n), nrow = m)), sparse = TRUE)
b <- rnorm(m)
I <- diag(n)
G <- rbind(cbind(I, -I),
           cbind(-I, -I))
G <- as(G, "dgCMatrix")
h <- rep(0.0, 2L * n)
dims <- list(l = 2L * n, q = NULL, e = 0L)
```

Note how ECOS expects sparse matrices, not ordinary matrices.

``` r

## Solve the problem
z <- ECOS_csolve(c = c, G = G, h = h, dims = dims, A = Atilde, b = b)
```

We check that the solution was found.

``` r

names(z)
```

    ## [1] "x"          "y"          "s"          "z"          "infostring"
    ## [6] "retcodes"   "summary"    "timing"

``` r

z$infostring
```

    ## [1] "Optimal solution found"

Extract the solution.

``` r

x <- z$x[1:n]
u <- z$x[(n+1):(2*n)]
nnzx = sum(abs(x) > 1e-8)
sprintf("x reconstructed with %d non-zero entries", nnzx / length(x) * 100)
```

    ## [1] "x reconstructed with 1 non-zero entries"

## Parametric Optimization with the Update-Solve Loop

When solving a sequence of related problems that share the same sparsity
structure but differ in numerical data, the multi-step lifecycle API
avoids repeating the expensive symbolic analysis and KKT ordering phase.
The pattern is:

1.  [`ECOS_setup()`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md)
    — one-time setup (symbolic analysis)
2.  Loop:
    [`ECOS_update()`](https://bnaras.github.io/ECOSolveR/reference/ECOS_update.md) +
    [`ECOS_solve()`](https://bnaras.github.io/ECOSolveR/reference/ECOS_solve.md)
    — cheap re-solves
3.  [`ECOS_cleanup()`](https://bnaras.github.io/ECOSolveR/reference/ECOS_cleanup.md)
    — free workspace

Here we vary the right-hand side $`h`$ of the LP from the previous
example, tracing how the optimal value changes as we scale $`h`$.

``` r

## Set up workspace once
ws <- ECOS_setup(c = c, G = G, h = h, dims = dims, A = Atilde, b = b)

## Sweep a parameter: scale h from 0 to 1
alphas <- seq(0, 1, length.out = 11)
pcosts <- numeric(length(alphas))

for (i in seq_along(alphas)) {
    ECOS_update(ws, h = h + alphas[i])
    res <- ECOS_solve(ws)
    pcosts[i] <- res$summary[["pcost"]]
}

ECOS_cleanup(ws)

## Show results
data.frame(alpha = alphas, pcost = round(pcosts, 4))
```

    ##    alpha     pcost
    ## 1    0.0   12.8872
    ## 2    0.1  -87.1128
    ## 3    0.2 -187.1128
    ## 4    0.3 -287.1128
    ## 5    0.4 -387.1128
    ## 6    0.5 -487.1128
    ## 7    0.6 -587.1128
    ## 8    0.7 -687.1128
    ## 9    0.8 -787.1128
    ## 10   0.9 -887.1128
    ## 11   1.0 -987.1128

For comparison, the same sweep using
[`ECOS_csolve()`](https://bnaras.github.io/ECOSolveR/reference/ECOS_csolve.md)
in a loop would call `ECOS_setup` and `ECOS_cleanup` at every iteration.
The lifecycle API is particularly beneficial when the problem dimensions
are large and setup dominates solve time.
