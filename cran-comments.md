## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

The CRAN incoming check flags "possibly misspelled words" in DESCRIPTION:
EVT and SAR are acronyms defined in the text (extreme value theory (EVT),
species-area relationship (SAR)); "et"/"al" are citation fragments; "endemism"
and "nestedness" are standard ecological terms. All are spelled correctly.

## Test environments

- Windows 11 (local), R 4.6.0
- win-builder: R-devel, R-release
- GitHub Actions: Ubuntu (R-release, R-devel), macOS (R-release), Windows (R-release)

## First submission

This is the first CRAN submission of spacc.

The package contains compiled C++ code (Rcpp/RcppParallel). The "New submission"
NOTE is expected. Examples and vignettes that depend on Suggests packages are
guarded with requireNamespace(). Parallel workers are capped to two cores when
_R_CHECK_LIMIT_CORES_ is set, so checks respect the CRAN core limit.
