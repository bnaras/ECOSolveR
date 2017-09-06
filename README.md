ECOSolveR
=========

Install this package the usual way in R using CRAN or via:

```{r}
devtools::install_github("bnaras/ECOSolveR")
```

This is an R wrapper around the [ecos](https://github.com/embotech/ecos) project
on GitHub.

Note that the ecos sources are included here. The only changes to the
ecos library source is to the makefiles; that was necessary to pass
the CRAN portability checks.

## Changes

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




