# Distance-Decay Analysis

Analyze how species richness changes with distance from focal points.

## Usage

``` r
distanceDecay(
  x,
  coords,
  n_seeds = 50L,
  breaks = NULL,
  distance = c("euclidean", "haversine"),
  progress = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  A site-by-species matrix.

- coords:

  A data.frame with x and y columns.

- n_seeds:

  Integer. Number of focal points.

- breaks:

  Numeric vector. Distance thresholds. Default auto-calculated.

- distance:

  Character. Distance method.

- progress:

  Logical. Show progress?

- seed:

  Integer. Random seed.

## Value

An object of class `spacc_decay` with distance-species relationships.

## References

Nekola, J.C. & White, P.S. (1999). The distance decay of similarity in
biogeography and ecology. Journal of Biogeography, 26, 867-878.

Soininen, J., McDonald, R. & Hillebrand, H. (2007). The distance decay
of similarity in ecological communities. Ecography, 30, 3-12.
