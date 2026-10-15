# Clean up an ECOS workspace

Frees the C-level ECOS workspace. After this call the workspace external
pointer is invalidated and cannot be used for further solves or updates.
It is safe to call this more than once.

## Usage

``` r
ECOS_cleanup(workspace)
```

## Arguments

- workspace:

  an external pointer of class `"ecos_workspace"` as returned by
  [`ECOS_setup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md).

## Value

`NULL`, invisibly.

## See also

[`ECOS_setup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md),
[`ECOS_solve`](https://bnaras.github.io/ECOSolveR/reference/ECOS_solve.md),
[`ECOS_update`](https://bnaras.github.io/ECOSolveR/reference/ECOS_update.md)
