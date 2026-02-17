# Reverse Dependency Check Results

**Date:** 2026-02-16
**Branch:** fix/audit-findings (commit 30732d7)
**R version:** 4.5.2 (2025-10-31)
**Platform:** aarch64-apple-darwin20 (macOS Tahoe 26.3)

## Summary

8 reverse dependencies checked. 0 failures attributable to ECOSolveR changes.

## Results

| Package | Version | Status | Notes |
|---------|---------|--------|-------|
| CVXR | 1.0-15 | 3 ERRORs | OSQP solver S7 `$` vs `@` incompatibility (confirmed: `solver="ECOS"` works fine) |
| EmpiricalDynamics | 0.1.2 | OK | |
| highOrderPortfolios | 0.1.1 | OK | |
| MOSAlloc | 1.2.5 | OK | |
| optiSel | 2.0.9 | 1 WARNING | `rgl.init` headless display warning |
| ROI.plugin.ecos | 1.0-2 | OK | |
| scpi | 3.0.1 | OK | |
| WpProj | 0.2.3 | 1 ERROR | Missing `bigmemory.sri` package |

## Details

### CVXR 1.0-15 (3 ERRORs -- not ECOSolveR-related)

The errors all stem from S7 0.2.1 breaking `solver$Solve` property access
in CVXR's OSQP solver class. When ECOS is used explicitly
(`solve(prob, solver="ECOS")`), CVXR works correctly. The default solver
selection picks OSQP for certain problem types, triggering the S7 error.
ECOSolveR is not mentioned anywhere in the error log.

### optiSel 2.0.9 (1 WARNING -- not ECOSolveR-related)

Warning from `rgl.init` failing on headless display. All tests and
examples pass.

### WpProj 0.2.3 (1 ERROR -- not ECOSolveR-related)

Missing `bigmemory.sri` package dependency. Tests pass; only one
example fails.

## Conclusion

The lifecycle API additions (`ECOS_setup`, `ECOS_solve`, `ECOS_update`,
`ECOS_cleanup`) are fully backward-compatible. All packages that directly
exercise ECOSolveR (especially ROI.plugin.ecos) pass without issues.
