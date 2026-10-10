## R CMD check results

0 errors | 0 warnings | 0 notes

This release adds richness estimators, coverage-standardised beta diversity,
network, redundancy and trait-variation accumulation, and moves every
accumulation method onto shared site orders in the C++ backend. The methods
are cited in the Description field.

It also fixes the NOTE on the current CRAN checks ("Rd files without
\usage: autoplot.spacc.Rd").

The bundled nanoflann header (inst/include/nanoflann.hpp, BSD license) has
been in the package since 0.8.3; its three copyright holders are now listed
in Authors@R as contributors and copyright holders, with the license text in
inst/COPYRIGHTS.

win-builder (R-devel and R-release) reports one NOTE listing possibly
misspelled words in DESCRIPTION: these are author surnames from the method
references, "EvoHeritage" (the name of the Rosindell et al. 2024 metric) and
"intraspecific", and are spelled correctly.

## Test environments

- Windows 11 (local), R 4.6.0
- win-builder: R-devel, R-release
- GitHub Actions: Ubuntu (R-release, R-devel), macOS (R-release), Windows (R-release)

## Downstream dependencies

No reverse dependencies.
