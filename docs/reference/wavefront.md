# Wavefront Expansion Accumulation

Accumulate species within an expanding radius from seed points. Models
invasion spread from introduction points.

## Usage

``` r
wavefront(
  x,
  coords,
  n_seeds = 50L,
  r0 = 0,
  dr = NULL,
  n_steps = 50L,
  distance = c("euclidean", "haversine"),
  progress = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  A site-by-species matrix.

- coords:

  A data.frame with x and y columns, or a `spacc_dist` object.

- n_seeds:

  Integer. Number of random starting points.

- r0:

  Numeric. Initial radius. Default 0.

- dr:

  Numeric. Radius increment per step. Default auto-calculated.

- n_steps:

  Integer. Number of expansion steps. Default 50.

- distance:

  Character. Distance method.

- progress:

  Logical. Show progress?

- seed:

  Integer. Random seed.

## Value

An object of class `spacc_wavefront` containing:

- curves:

  Matrix of species counts (n_seeds x n_steps)

- radius:

  Vector of radius values

- sites_included:

  Matrix of sites included at each step

## Examples

``` r
if (FALSE) { # \dontrun{
wf <- wavefront(species, coords, n_seeds = 20, n_steps = 100)
plot(wf)
} # }
```
