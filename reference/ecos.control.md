# Return the default optimization parameters for ECOS

This is used to control the behavior of the underlying optimization
code.

## Usage

``` r
ecos.control(
  maxit = 100L,
  feastol = 1e-08,
  reltol = 1e-08,
  abstol = 1e-08,
  feastol_inacc = 1e-04,
  abstol_inacc = 5e-05,
  reltol_inacc = 5e-05,
  verbose = 0L,
  mi_max_iters = 1000L,
  mi_int_tol = 1e-04,
  mi_abs_eps = 1e-06,
  mi_rel_eps = 1e-06
)
```

## Arguments

- maxit:

  the maximum number of iterations for ecos, default 100L

- feastol:

  the tolerance on the primal and dual residual, default 1e-8

- reltol:

  the relative tolerance on the duality gap, default 1e-8

- abstol:

  the absolute tolerance on the duality gap, default 1e-8

- feastol_inacc:

  the tolerance on the primal and dual residual if reduced precisions,
  default 1e-4

- abstol_inacc:

  the absolute tolerance on the duality gap if reduced precision,
  default 5e-5

- reltol_inacc:

  the relative tolerance on the duality gap if reduced precision,
  default 5e-5

- verbose:

  verbosity level, default 0L. A verbosity level of 1L will show more
  detail, but clutter session transcript.

- mi_max_iters:

  the maximum number of branch and bound iterations (mixed integer
  problems only), default 1000L

- mi_int_tol:

  the integer tolerence (mixed integer problems only), default 1e-4

- mi_abs_eps:

  the absolute tolerance between upper and lower bounds (mixed integer
  problems only), default 1e-6

- mi_rel_eps:

  the relative tolerance, \\(U-L)/L\\, between upper and lower bounds
  (mixed integer problems only), default 1e-6

## Value

a list with the following elements:

- FEASTOL:

  the tolerance on the primal and dual residual, parameter `feastol`

- ABSTOL:

  the absolute tolerance on the duality gap, parameter `abstol`

- RELTOL:

  the relative tolerance on the duality gap, parameter `reltol`

- FEASTOL_INACC:

  the tolerance on the primal and dual residual if reduced precisions,
  parameter `feastol_inacc`

- ABSTOL_INACC:

  the absolute tolerance on the duality gap if reduced precision,
  parameter `abstol_inacc`

- RELTOL_INACC:

  the relative tolerance on the duality gap if reduced precision,
  parameter `reltol_inacc`

- MAXIT:

  the maximum number of iterations for ecos, parameter `maxit`

- MI_MAX_ITERS:

  the maximum number of branch and bound iterations (mixed integer
  problems only), parameter `mi_max_iters`

- MI_INT_TOL:

  the integer tolerence (mixed integer problems only), parameter
  `mi_int_tol`

- MI_ABS_EPS:

  the absolute tolerance between upper and lower bounds (mixed integer
  problems only), parameter `mi_abs_eps`

- MI_REL_EPS:

  the relative tolerance, \\(U-L)/L\\, between upper and lower bounds
  (mixed integer problems only), parameter `mi_rel_eps`

- VERBOSE:

  verbosity level, parameter `verbose`
