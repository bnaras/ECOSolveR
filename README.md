
<!-- README.md is generated from README.Rmd. Please edit that file -->


# ECOSolveR

  [![R-CMD-check](https://github.com/bnaras/ECOSolveR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bnaras/ECOSolveR/actions/workflows/R-CMD-check.yaml)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/ECOSolveR)](https://cran.r-project.org/package=ECOSolveR)
[![Coverage
Status](https://img.shields.io/codecov/c/github/bnaras/ECOSolveR/master.svg)](https://app.codecov.io/github/bnaras/ECOSolveR?branch=master)
[![](https://cranlogs.r-pkg.org/badges/ECOSolveR)](https://CRAN.R-project.org/package=ECOSolveR)

Embedded Conic Solver in R. This is an R wrapper around the
[ecos](https://github.com/embotech/ecos) project on GitHub which
describes ECOS as below.

ECOS is a numerical software for solving convex second-order cone
programs (SOCPs) of type

$$
\mbox{Minimize } c'x \mbox{ such that } {\mathbf Ax} = {\mathbf b} \mbox{ and } {\mathbf G \mathbf x}\,\, \leq_{\mathbf K}\,\, {\mathbf h}
$$
where the last inequality is generalized, that is, ${\mathbf h}-\mathbf{Gx}$ belongs to
the cone ${\mathbf K}$.

ECOS supports the positive orthant ${\mathbf R}_+$, second-order cones
${\mathbf Q}_n$ defined as

$$
{\mathbf Q}_n = \bigl\{ (t,{\mathbf x}) | t >= \lVert{\mathbf x}\rVert_2 \bigr\}
$$

with $t$ a scalar and ${\mathbf x} \in {\mathbf R}_{n-1}$, and the exponential
cone ${\mathbf K}_e$ defined as

$$
\mathbf{K}_e = \mbox{closure} \bigl\{ (x,y,z) | exp(x/z) <= y/z, z>0 \bigr\},
$$

where  $(x,y,z) \in {\mathbf R}^3$.

The cone ${\mathbf K}$ is therefore a direct product of the positive orthant, second-order, and exponential cones:

$$
{\mathbf K} = {\mathbf R}_+ \times {\mathbf Q}_{n_1} \times \cdots \times {\mathbf Q}_{n_N} \times {\mathbf K}_e \times \cdots \times {\mathbf K}_e.
$$

## Further Details

Note that the ECOS C language sources are included here. Changes to
the original source are clearly delineated for easy reference.



