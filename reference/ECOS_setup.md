# Set up an ECOS workspace for multi-step solving

Creates an ECOS workspace that can be solved, updated with new numerical
data, and solved again without repeating the expensive symbolic analysis
phase. This is useful for parametric optimization, model predictive
control, and other settings where the problem structure stays the same
but data changes.

## Usage

``` r
ECOS_setup(
  c,
  G,
  h,
  dims = list(l = integer(0), q = NULL, e = integer(0)),
  A = NULL,
  b = numeric(0),
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

- control:

  is a named list that controls various optimization parameters; see
  [ecos.control](https://bnaras.github.io/ECOSolveR/reference/ecos.control.md).

## Value

an external pointer of class `"ecos_workspace"`. Must eventually be
freed via
[`ECOS_cleanup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_cleanup.md)
or R garbage collection.

## See also

[`ECOS_solve`](https://bnaras.github.io/ECOSolveR/reference/ECOS_solve.md),
[`ECOS_update`](https://bnaras.github.io/ECOSolveR/reference/ECOS_update.md),
[`ECOS_cleanup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_cleanup.md),
[`ECOS_csolve`](https://bnaras.github.io/ECOSolveR/reference/ECOS_csolve.md)
