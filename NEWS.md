# ECOSolveR 0.5-3

* Address the issue caused by Matrix 1.3.x: Force coercion to
  "dgCMatrix" if not so. Make corresponding change in vignette.

# ECOSolveR 0.5-3

* Fix some awful cut-and-paste error in `R/ecos.R`
* Address edge case for ECOS_BB when `G` is `NULL`.

# ECOSolveR 0.5-2

* Fix compilation issue on Ubuntu by modifying `Makevars` and
  `Makevars.win` to use `Rscript` call to `R.home()` to figure out
  include path. ([Issue 5](https://github.com/bnaras/ECOSolveR/issues/5))

# ECOSolveR 0.5-1

* ECOS source fix for header `glblopts.h` that defined `ECOS_NAN`
  using portable R version `R_NaN` and `R_PosInf`. This caused the
  convolution example in `CVXR` package to fail on 32-bit
  platforms. 

* Added a convolution test to be specific.

# ECOSolveR 0.5

* ECOS update: Synced underlying library to version 2.0.7.

* Removed import of `Matrix` package, added `slam` interface,
  contributions from Florian Schwendinger.

* Tests: Added a number of unit tests based on base C library.

# ECOSolveR 0.4

* `ECOS_csolve` assumes `A` and `G` as `NULL` if any dimension is 0
* The dims list is now coerced to integer immediately upon entry to
  `ECOS_csolve`
* ECOS library sources updated to version 2.0.5.

# ECOSolveR 0.3-[1-3]

* Updated Readme
* Allowed empty vector `c`, per Anqi’s request

# ECOSolveR 0.3

* Registered `.Call` entries
* Changed Anqi’s email

# ECOSolveR 0.2

* Added system requirement of GNU makefile
* Improved the description of the package
* Checked more carefully on Linux, MacOS, and Windows

# ECOSolveR 0.1-1

* First version, thought I had it right for all platforms, but quite
  wrong. 

