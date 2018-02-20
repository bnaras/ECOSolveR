ECOSolveR
=========

Embedded Conic Solver in R. This is an R wrapper around the
[ecos](https://github.com/embotech/ecos) project on GitHub which
describes ECOS as below.

ECOS is a numerical software for solving convex second-order cone programs (SOCPs) of type

\begin{array}{llll}
\mbox{Minimize} & c'x & \ni Ax = b & \mbox{and} & Gx <=_K h
\end{array}

where the last inequality is generalized, i.e. $ h - Gx $ belongs to the
cone $ K. $ 

ECOS supports the positive orthant $ {\mathbf R}_+ $, second-order cones
$ Q_n $ defined as

$$
Q_n = { (t,x) | t >= || x ||_2 }
$$

with $t$ a scalar and $x \in {\mathbf R}_{n-1}$, and the exponential
cone $K_e$ defined as

$$
K_e = {\mbox{closure}} \{ (x,y,z) | exp(x/z) <= y/z, z>0 \},
$$

where  $(x,y,z) \in {\mathbf R}^3$.

The cone $K$ is therefore a direct product of the positive orthant, second-order, and exponential cones:

\begin{array}{l}
K = R_+ \times Q_{n_1} \times \cdots \times Q_{n_N} \times K_e \times \cdots \times K_e.
\end{array}

## Further Details

Note that the ECOS C language sources are included here. The only
edits to the ECOS library source is to the makefiles; that was
necessary to pass the CRAN portability checks.


## Changes

- Version 0.3-3 and 0.4 (2018-02-16)
	- ECOS_csolve assumes A and G as NULL if any dimension is 0
	- The dims list is now coerced to integer immediately upon
	  entry to ECOS_csolve
	- ECOS library sources updated to version 2.0.5.

- Version 0.3-1 and 0.3-2 (2017-09-05)
	- Updated Readme
	- Allowed empty vector c, per Anqi's request

- Version 0.3 (2017-05-16)
	- Registered .Call entries
	- Changed Anqi's email

- Version 0.2 (2015-11-17)
	- Added system requirement of GNU makefile
	- Improved the description of the package
	- Checked more carefully on Linux, MacOS, and Windows

- Version 0.1-1 (2015-11-17)
	- First version, thought I had it right for all platforms, but quite wrong!


