# Solve a conic optimization problem

The function `ECOS_csolve` is a wrapper around the ecos `csolve` C
function. Conic constraints are specified using the \\G\\ and \\h\\
parameters and can be `NULL` and zero length vector respectively
indicating an absence of conic constraints. Similarly, equality
constraints are specified via \\A\\ and \\b\\ parameters with `NULL` and
empty vector values representing a lack of such constraints. At most one
of the pair \\(G , h)\\ or \\(A, b)\\ is allowed to be absent.

## Usage

``` r
ECOS_csolve(
  c = numeric(0),
  G = NULL,
  h = numeric(0),
  dims = list(l = integer(0), q = NULL, e = integer(0)),
  A = NULL,
  b = numeric(0),
  bool_vars = integer(0),
  int_vars = integer(0),
  control = ecos.control()
)
```

## Arguments

- c:

  the coefficients of the objective function; the length of this
  determines the number of variables \\n\\ in the problem.

- G:

  the inequality constraint matrix in one of three forms: a plain
  matrix, simple triplet matrix, or compressed column format, e.g.
  [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).
  Can also be `NULL`

- h:

  the right hand side of the inequality constraint. Can be empty numeric
  vector.

- dims:

  is a list of three named elements: `dims['l']` an integer specifying
  the dimension of positive orthant cone, `dims['q']` an integer vector
  specifying dimensions of second-order cones, `dims['e']` an integer
  specifying the number of exponential cones

- A:

  the optional equality constraint matrix in one of three forms: a plain
  matrix, simple triplet matrix, or compressed column format, e.g.
  [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).
  Can be `NULL`

- b:

  the right hand side of the equality constraint, must be specified if
  \\A\\ is. Can be empty numeric vector.

- bool_vars:

  the indices of the variables, 1 through \\n\\, that are boolean; that
  is, they are either present or absent in the solution

- int_vars:

  the indices of the variables, 1 through \\n\\, that are integers

- control:

  is a named list that controls various optimization parameters; see
  [ecos.control](https://bnaras.github.io/ECOSolveR/reference/ecos.control.md).

## Value

a list of 8 named items

- x:

  primal variables

- y:

  dual variables for equality constraints

- s:

  slacks for \\Gx + s \<= h\\, \\s \in K\\

- z:

  dual variables for inequality constraints \\s \in K\\

- infostring:

  gives information about the status of solution

- retcodes:

  a named integer vector containing four elements

  exitflag

  :   0=`ECOS_OPTIMAL`, 1=`ECOS_PINF`, 2=`ECOS_DINF`,
      10=`ECOS_INACC_OFFSET`, -1=`ECOS_MAXIT`, -2=`ECOS_NUMERICS`,
      -3=`ECOS_OUTCONE`, -4=`ECOS_SIGINT`, -7=`ECOS_FATAL`. See
      [ECOS_exitcodes](https://bnaras.github.io/ECOSolveR/reference/ECOS_exitcodes.md)

  iter

  :   the number of iterations used

  mi_iter

  :   the number of iterations for mixed integer problems

  numerr

  :   a non-zero number if a numeric error occurred

- summary:

  a named numeric vector containing

  pcost

  :   value of primal objective

  dcost

  :   value of dual objective

  pres

  :   primal residual on inequalities and equalities

  dres

  :   dual residual

  pinf

  :   primal infeasibility measure

  dinf

  :   dual infeasibility measure

  pinfres

  :   primal infeasibility residual

  dinfres

  :   dual infeasibility residual

  gap

  :   duality gap

  relgap

  :   relative duality gap

  r0

  :   Unknown at the moment to this R package maintainer.

- timing:

  a named numeric vector of timing information consisting of

  runtime

  :   the total runtime in ecos

  tsetup

  :   the time for setup of the problem

  tsolve

  :   the time to solve the problem

## Details

A call to this function will solve the problem: minimize \\c^Tx\\,
subject to \\Ax = b\\, and \\h - G\*x \in K\\.

Variables can be constrained to be boolean (1 or 0) or integers. This is
indicated by specifying parameters `bool_vars` and/or `int_vars`
respectively. If so indicated, the solutions will be found using a
branch and bound algorithm.

## Examples

``` r

## githubIssue98
cat("Basic matrix interface\n")
#> Basic matrix interface
Gmat <- matrix(c(0.416757847405471, 2.13619609566845, 1.79343558519486, 0, 0,
                 0, 0, -1, 0, 0, 0, 0.056266827226329, -1.64027080840499, 0.841747365656204,
                 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0.416757847405471, 2.13619609566845,
                 1.79343558519486, 0, 0, 0, -1, 0, 0, 0, 0, 0.056266827226329, -1.64027080840499,
                 0.841747365656204, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0), ncol = 5L)
c <- as.numeric(c(0, 0, 0, 0, 1))
h <- as.numeric(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
dims <- list(l = 6L, q = 5L, e = 0L)
ECOS_csolve(c = c, G = Gmat, h = h,
           dims = dims,
           A = NULL, b = numeric(0))
#> $x
#> [1]  1.479062e-11  4.332909e-12  1.478861e-11  4.273179e-12 -3.372330e-10
#> 
#> $y
#> numeric(0)
#> 
#> $s
#>  [1] 5.007707e-11 2.112317e-10 1.720096e-10 5.008123e-11 2.111387e-10
#>  [6] 1.720630e-10 2.502815e-10 1.887750e-11 4.007416e-12 1.887598e-11
#> [11] 4.038335e-12
#> 
#> $z
#>  [1]  0.36678432  0.07561476  0.09859791  0.36678392  0.07561520  0.09859751
#>  [7]  1.00000000  0.49121720 -0.02039636  0.49121725 -0.02039743
#> 
#> $infostring
#> [1] "Optimal solution found"
#> 
#> $retcodes
#> exitFlag     iter  mi_iter   numerr 
#>        0        6       -1        0 
#> 
#> $summary
#>         pcost         dcost          pres          dres          pinf 
#> -3.372330e-10  0.000000e+00  4.264405e-10  5.482329e-11  0.000000e+00 
#>          dinf       pinfres       dinfres           gap        relgap 
#>  0.000000e+00           NaN           NaN  1.738783e-09  5.156029e+00 
#>            r0 
#>  1.000000e-08 
#> 
#> $timing
#>     runtime      tsetup      tsolve 
#> 0.000129975 0.000025588 0.000104387 
#> 

cat("Simple Triplet Matrix interface, if you have package slam\n")
#> Simple Triplet Matrix interface, if you have package slam
if (requireNamespace("slam")) {

  ECOS_csolve(c = c, G = slam::as.simple_triplet_matrix(Gmat), h = h,
              dims = dims,
              A = NULL, b = numeric(0))
}
#> Loading required namespace: slam
#> $x
#> [1]  1.479062e-11  4.332909e-12  1.478861e-11  4.273179e-12 -3.372330e-10
#> 
#> $y
#> numeric(0)
#> 
#> $s
#>  [1] 5.007707e-11 2.112317e-10 1.720096e-10 5.008123e-11 2.111387e-10
#>  [6] 1.720630e-10 2.502815e-10 1.887750e-11 4.007416e-12 1.887598e-11
#> [11] 4.038335e-12
#> 
#> $z
#>  [1]  0.36678432  0.07561476  0.09859791  0.36678392  0.07561520  0.09859751
#>  [7]  1.00000000  0.49121720 -0.02039636  0.49121725 -0.02039743
#> 
#> $infostring
#> [1] "Optimal solution found"
#> 
#> $retcodes
#> exitFlag     iter  mi_iter   numerr 
#>        0        6       -1        0 
#> 
#> $summary
#>         pcost         dcost          pres          dres          pinf 
#> -3.372330e-10  0.000000e+00  4.264405e-10  5.482329e-11  0.000000e+00 
#>          dinf       pinfres       dinfres           gap        relgap 
#>  0.000000e+00           NaN           NaN  1.738783e-09  5.156029e+00 
#>            r0 
#>  1.000000e-08 
#> 
#> $timing
#>    runtime     tsetup     tsolve 
#> 4.3044e-05 1.6815e-05 2.6229e-05 
#> 

if (requireNamespace("Matrix")) {
   ECOS_csolve(c = c, G = Matrix::Matrix(Gmat), h = h,
               dims = dims,
               A = NULL, b = numeric(0))
}
#> Loading required namespace: Matrix
#> $x
#> [1]  1.479062e-11  4.332909e-12  1.478861e-11  4.273179e-12 -3.372330e-10
#> 
#> $y
#> numeric(0)
#> 
#> $s
#>  [1] 5.007707e-11 2.112317e-10 1.720096e-10 5.008123e-11 2.111387e-10
#>  [6] 1.720630e-10 2.502815e-10 1.887750e-11 4.007416e-12 1.887598e-11
#> [11] 4.038335e-12
#> 
#> $z
#>  [1]  0.36678432  0.07561476  0.09859791  0.36678392  0.07561520  0.09859751
#>  [7]  1.00000000  0.49121720 -0.02039636  0.49121725 -0.02039743
#> 
#> $infostring
#> [1] "Optimal solution found"
#> 
#> $retcodes
#> exitFlag     iter  mi_iter   numerr 
#>        0        6       -1        0 
#> 
#> $summary
#>         pcost         dcost          pres          dres          pinf 
#> -3.372330e-10  0.000000e+00  4.264405e-10  5.482329e-11  0.000000e+00 
#>          dinf       pinfres       dinfres           gap        relgap 
#>  0.000000e+00           NaN           NaN  1.738783e-09  5.156029e+00 
#>            r0 
#>  1.000000e-08 
#> 
#> $timing
#>    runtime     tsetup     tsolve 
#> 5.7027e-05 2.9064e-05 2.7963e-05 
#> 

## Larger problems using saved data can be found in the test suite.
## Here is one
if (requireNamespace("Matrix")) {
  MPC01 <- readRDS(system.file("testdata", "MPC01_1.RDS", package = "ECOSolveR"))
  G <- Matrix::sparseMatrix(x = MPC01$Gpr, i = MPC01$Gir, p = MPC01$Gjc,
                            dims = c(MPC01$m, MPC01$n), index1 = FALSE)
  h <- MPC01$h
  dims <- lapply(list(l = MPC01$l, q=MPC01$q, e=MPC01$e), as.integer)
  retval <- ECOS_csolve(c = MPC01$c, G=G, h = h, dims = dims, A = NULL, b = NULL,
                        control = ecos.control(verbose=1L))
  retval$retcodes
  retval$infostring
  retval$summary
}
#> 
#> ECOS 2.0.7 - (C) embotech GmbH, Zurich Switzerland, 2012-15. Web: www.embotech.com/ECOS
#> 
#> It     pcost       dcost      gap   pres   dres    k/t    mu     step   sigma     IR    |   BT
#>  0  -1.030e+02  -3.786e+03  +6e+04  8e-01  3e-01  1e+00  3e+01    ---    ---    1  1  - |  -  - 
#>  1  -1.873e+02  -5.983e+02  +1e+04  7e-02  3e-02  5e-01  6e+00  0.8442  3e-02   1  2  2 |  0  0
#>  2  -1.890e+02  -2.683e+02  +4e+03  1e-02  6e-03  4e-01  2e+00  0.8197  2e-01   2  3  3 |  0  0
#>  3  -2.094e+02  -2.160e+02  +4e+02  1e-03  4e-04  3e-02  2e-01  0.9129  2e-02   2  2  2 |  0  0
#>  4  -2.109e+02  -2.113e+02  +3e+01  7e-05  2e-05  6e-04  2e-02  0.9263  3e-03   2  3  3 |  0  0
#>  5  -2.112e+02  -2.112e+02  +8e-01  1e-06  5e-07  1e-05  4e-04  0.9779  4e-03   4  7  8 |  0  0
#>  6  -2.112e+02  -2.112e+02  +9e-03  2e-08  6e-09  2e-07  5e-06  0.9890  1e-04   2  1  1 |  0  0
#>  7  -2.112e+02  -2.112e+02  +1e-04  2e-10  6e-11  2e-09  6e-08  0.9890  1e-04   2  2  2 |  0  0
#>  8  -2.112e+02  -2.112e+02  +1e-06  2e-12  7e-13  2e-11  6e-10  0.9890  1e-04   2  1  2 |  0  0
#> 
#> OPTIMAL (within feastol=2.4e-12, reltol=5.3e-09, abstol=1.1e-06).
#> Runtime: 0.127031 seconds.
#> 
#>         pcost         dcost          pres          dres          pinf 
#> -2.112370e+02 -2.112370e+02  2.387958e-12  7.134616e-13  0.000000e+00 
#>          dinf       pinfres       dinfres           gap        relgap 
#>  0.000000e+00           NaN  4.964435e-01  1.111838e-06  5.263464e-09 
#>            r0 
#>  1.000000e-08 
```
