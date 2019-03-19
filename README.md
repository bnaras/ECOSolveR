
<!-- README.md is generated from README.Rmd. Please edit that file -->



# ECOSolveR

[![Travis-CI Build Status](https://travis-ci.org/bnaras/ECOSolveR.svg?branch=master)](https://travis-ci.org/bnaras/ECOSolveR)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/ECOSolveR)](https://cran.r-project.org/package=ECOSolveR)
[![Coverage Status](https://codecov.io/gh/bnaras/ECOSolveR/branch/master/graph/badge.svg)](https://codecov.io/gh/bnaras/ECOSolveR)

Embedded Conic Solver in R. This is an R wrapper around the
[ecos](https://github.com/embotech/ecos) project on GitHub which
describes ECOS as below.

ECOS is a numerical software for solving convex second-order cone programs (SOCPs) of type

$$
\mbox{Minimize } c'x \mbox{ such that } Ax = b \mbox{ and } Gx <=_K h
$$
where the last inequality is generalized, that is, $h-Gx$ belongs to
the cone $K$.

ECOS supports the positive orthant ${\mathbf R}_+$, second-order cones
$Q_n$ defined as

$$
Q_n = \bigl\{ (t,{\mathbf x}) | t >= \lVert{\mathbf x}\rVert_2 \bigr\}
$$

with $t$ a scalar and ${\mathbf x} \in {\mathbf R}_{n-1}$, and the exponential
cone $K_e$ defined as

$$
K_e = \mbox{closure} \bigl\{ (x,y,z) | exp(x/z) <= y/z, z>0 \bigr\},
$$

where  $(x,y,z) \in {\mathbf R}^3$.

The cone $K$ is therefore a direct product of the positive orthant, second-order, and exponential cones:

$$
\begin{array}{l}
K = R_+ \times Q_{n_1} \times \cdots \times Q_{n_N} \times K_e \times \cdots \times K_e.
\end{array}
$$

## Further Details

Note that the ECOS C language sources are included here. Changes to
the original source are clearly delineated for easy reference.



