# Contributing to spacc

## Code of conduct

Participation in this project follows
[CODE_OF_CONDUCT.md](https://gillescolling.com/spacc/CODE_OF_CONDUCT.md).

## Setup

Clone the repository and install the development dependencies from R:

``` r

pak::pak(c("devtools", "testthat", "rcmdcheck", "roxygen2"))
devtools::install_deps(dependencies = TRUE)
```

## Testing

Run the test suite and package check before opening a pull request:

``` r

devtools::test()
devtools::check()
```

Changes to statistical methods need tests against independently
calculated reference values. Changes to estimators need
parameter-recovery or interval- coverage tests using simulated data with
known truth.

## Documentation

Document exported functions with roxygen2 and regenerate the package
files:

``` r

devtools::document()
```

Update the relevant vignette and `NEWS.md` when behavior visible to
users changes.

## Project organization

- `R/` contains the R interface and S3 methods.
- `src/` contains the compiled accumulation kernels.
- `tests/testthat/` contains unit, reference, recovery, and coverage
  tests.
- `vignettes/` contains worked documentation.
- `paper/` contains the manuscript and its figures.

## Pull requests

1.  Create a focused branch.
2.  Add the implementation and tests in the same change.
3.  Run
    [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
    and
    [`devtools::check()`](https://devtools.r-lib.org/reference/check.html).
4.  Regenerate documentation when public interfaces change.
5.  Describe the behavior change and the evidence supporting it.

Contributions are licensed under the project license when merged.
