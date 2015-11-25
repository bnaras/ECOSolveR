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