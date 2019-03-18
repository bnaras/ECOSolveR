#' ECOS solver exit codes
#'
#' A two-column data frame consisting of the code and description for the
#' ECOS solver with ECOS symbolic code names as row names
#'
#' @name ECOS_exitcodes
#' @docType data
#' @keywords data
#'
NULL


## Source: ecos/include/ecos.h
## ECOS_exitcodes <- data.frame(
##     code = c(0L, 1L, 2L, 10L, -1L, -2L, -3L, -4L, -7L),
##     description = c("Problem solved to optimality",
##                     "Found certificate of primal infeasibility",
##                     "Found certificate of dual infeasibility",
##                     "Offset exitflag at inaccurate results",
##                     "Maximum number of iterations reached",
##                     "Search direction unreliable",
##                     "s or z got outside the cone, numerics?",
##                     "solver interrupted by a signal/ctrl-c",
##                     "Unknown problem in solver"),
##     stringsAsFactors = FALSE)
## rownames(ECOS_exitcode) <-
##             c("ECOS_OPTIMAL",
##              "ECOS_PINF",
##              "ECOS_DINF",
##              "ECOS_INACC_OFFSET",
##              "ECOS_MAXIT",
##              "ECOS_NUMERICS",
##              "ECOS_OUTCONE",
##              "ECOS_SIGINT",
##              "ECOS_FATAL")
##
##


