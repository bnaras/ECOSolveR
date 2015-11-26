ECOSolveR
=========

Install this package the usual way in R once it gets on CRAN (very shortly) or via:

```{r}
library(devtools)
install_github("bnaras/ECOSolveR")
```

This is an R wrapper around the [ecos](https://github.com/embotech/ecos) project
on GitHub.

Note that the ecos sources are included here. The only changes I have made are to
the ecos library stuff is to the makefiles; this was necessary to pass the
CRAN portability checks.

## Changes

- Version 0.1-1 (2015-11-17)
	- First version, thought I had it right for all platforms, but quite wrong!

- Version 0.2 (2015-11-17)
	- Added system requirement of GNU makefile
	- Improved the description of the package
	- Checked more carefully on Linux, MacOS, and Windows
