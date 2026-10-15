# Update numerical data in an ECOS workspace

Replaces numerical values in an existing ECOS workspace without
repeating symbolic analysis. The sparsity structure must remain the
same; only the non-zero values of G and A, and the dense vectors c, h, b
can be updated. Pass `NULL` for any argument to leave it unchanged.

## Usage

``` r
ECOS_update(workspace, Gpr = NULL, Apr = NULL, c = NULL, h = NULL, b = NULL)
```

## Arguments

- workspace:

  an external pointer of class `"ecos_workspace"` as returned by
  [`ECOS_setup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md).

- Gpr:

  numeric vector of new non-zero values for G (same length as original),
  or `NULL` to keep current values.

- Apr:

  numeric vector of new non-zero values for A (same length as original),
  or `NULL` to keep current values.

- c:

  numeric vector of new objective coefficients, or `NULL`.

- h:

  numeric vector of new inequality RHS, or `NULL`.

- b:

  numeric vector of new equality RHS, or `NULL`.

## Value

the `workspace` object, invisibly.

## See also

[`ECOS_setup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md),
[`ECOS_solve`](https://bnaras.github.io/ECOSolveR/reference/ECOS_solve.md),
[`ECOS_cleanup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_cleanup.md)
