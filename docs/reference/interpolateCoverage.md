# Interpolate Richness at Target Coverage Levels

Estimate species richness at specified coverage levels by interpolation.

## Usage

``` r
interpolateCoverage(x, target = c(0.9, 0.95, 0.99))
```

## Arguments

- x:

  A `spacc_coverage` object.

- target:

  Numeric vector of target coverage levels (0 to 1). Default
  `c(0.90, 0.95, 0.99)`.

## Value

A data.frame with columns for each target coverage level, showing
interpolated richness for each seed.
