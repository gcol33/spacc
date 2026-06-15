# Predict richness at a given sampling effort

Predict richness at a given sampling effort

## Usage

``` r
# S3 method for class 'spacc_fit'
predict(object, n = NULL, interval = c("none", "bootstrap"), level = NULL, ...)
```

## Arguments

- object:

  A `spacc_fit` object.

- n:

  Integer vector of site counts. Defaults to the observed range.

- interval:

  `"none"` (default) returns point predictions; `"bootstrap"` returns a
  data frame with percentile bounds from the stored bootstrap draws.

- level:

  Confidence level for the interval. Defaults to the fit's level.

- ...:

  Unused.

## Value

A numeric vector (point) or a data frame with `n`, `fit`, `lower`,
`upper` (interval).
