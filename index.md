# ECOSolveR

[![R-CMD-check](https://github.com/bnaras/ECOSolveR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bnaras/ECOSolveR/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/ECOSolveR)](https://cran.r-project.org/package=ECOSolveR)
[![Coverage
Status](https://img.shields.io/codecov/c/github/bnaras/ECOSolveR/master.svg)](https://app.codecov.io/github/bnaras/ECOSolveR?branch=master)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/ECOSolveR)](https://CRAN.R-project.org/package=ECOSolveR)

Embedded Conic Solver in R. This is an R wrapper around the
[ecos](https://github.com/embotech/ecos) project on GitHub which
describes ECOS as below.

ECOS is a numerical software for solving convex second-order cone
programs (SOCPs) of type

\text{Minimize } \mathbf{c}^\top\mathbf{x} \text{ such that }
\mathbf{A}\mathbf{x} = \mathbf{b} \text{ and } \mathbf{G}\mathbf{x}
\leq\_{\mathbf{K}} \mathbf{h}

where the last inequality is generalized, that is, \mathbf{h} -
\mathbf{G}\mathbf{x} belongs to the cone \mathbf{K}.

ECOS supports the positive orthant \mathbf{R}\_+, second-order cones
\mathbf{Q}\_n defined as

\mathbf{Q}\_n = \bigl\lbrace (t, \mathbf{x}) \mid t \geq \lVert
\mathbf{x} \rVert_2 \bigr\rbrace

with t a scalar and \mathbf{x} \in \mathbf{R}^{n-1}, and the exponential
cone \mathbf{K}\_e defined as

\mathbf{K}\_e = \text{closure} \bigl\lbrace (x, y, z) \mid \exp(x/z)
\leq y/z,\\ z \> 0 \bigr\rbrace,

where (x, y, z) \in \mathbf{R}^3.

The cone \mathbf{K} is therefore a direct product of the positive
orthant, second-order, and exponential cones:

\mathbf{K} = \mathbf{R}\_+ \times \mathbf{Q}\_{n_1} \times \cdots \times
\mathbf{Q}\_{n_N} \times \mathbf{K}\_e \times \cdots \times
\mathbf{K}\_e.

## Further Details

Note that the ECOS C language sources are included here. Changes to the
original source are clearly delineated for easy reference.
