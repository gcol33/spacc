# Predict richness at a given area

Predict richness at a given area

## Usage

``` r
# S3 method for class 'spacc_area'
predict(
  object,
  area = NULL,
  interval = c("none", "bootstrap"),
  level = NULL,
  ...
)
```

## Arguments

- object:

  A `spacc_area` object.

- area:

  Numeric vector of areas. Defaults to the observed T-S areas.

- interval:

  `"none"` (default) or `"bootstrap"` for percentile bounds.

- level:

  Confidence level. Defaults to the fit's level.

- ...:

  Unused.

## Value

A numeric vector (point) or a data frame with `area`, `fit`, `lower`,
`upper`.
