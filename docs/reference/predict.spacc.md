# Predict from Empirical Accumulation Curve

Interpolate the mean empirical accumulation curve at arbitrary site
counts using linear interpolation. Unlike the predict method for
`spacc_fit` objects, this does not use a fitted model; it interpolates
the observed curve directly.

## Usage

``` r
# S3 method for class 'spacc'
predict(object, n = NULL, ci = TRUE, ci_level = 0.95, warn = TRUE, ...)
```

## Arguments

- object:

  A `spacc` object.

- n:

  Numeric vector of site counts at which to interpolate. Defaults to
  25%, 50%, and 100% of total sites.

- ci:

  Logical. If `TRUE` (default), return a data frame with columns `n`,
  `mean`, `lower`, `upper`. If `FALSE`, return a named numeric vector.

- ci_level:

  Confidence level for the interval (default 0.95).

- warn:

  Logical. Warn when `n` values fall outside the observed range (default
  `TRUE`).

- ...:

  Ignored.

## Value

A data frame (if `ci = TRUE`) or named numeric vector (if `ci = FALSE`).
Out-of-range values return `NA`.

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
sac <- spacc(species, coords, n_seeds = 10, progress = FALSE)
predict(sac, n = c(10, 25, 50))
predict(sac, n = c(10, 25), ci = FALSE)
# }
```
