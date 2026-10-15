# Solve an ECOS workspace

Calls `ECOS_solve` on an existing workspace created by
[`ECOS_setup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md).
The workspace can be solved multiple times, optionally with
[`ECOS_update`](https://bnaras.github.io/ECOSolveR/reference/ECOS_update.md)
calls in between to change numerical data.

## Usage

``` r
ECOS_solve(workspace, control = NULL)
```

## Arguments

- workspace:

  an external pointer of class `"ecos_workspace"` as returned by
  [`ECOS_setup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md).

- control:

  optional named list of solver parameters (see
  [`ecos.control`](https://bnaras.github.io/ECOSolveR/reference/ecos.control.md)).
  If `NULL` (the default), the settings from
  [`ECOS_setup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md)
  (or the most recent `ECOS_solve` call with non-NULL control) are
  reused.

## Value

the same result list as
[`ECOS_csolve`](https://bnaras.github.io/ECOSolveR/reference/ECOS_csolve.md).

## See also

[`ECOS_setup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_setup.md),
[`ECOS_update`](https://bnaras.github.io/ECOSolveR/reference/ECOS_update.md),
[`ECOS_cleanup`](https://bnaras.github.io/ECOSolveR/reference/ECOS_cleanup.md)
