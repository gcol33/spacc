# Spatial Subsampling

Reduce spatial autocorrelation by subsampling sites using various
methods.

## Usage

``` r
subsample(
  coords,
  n = NULL,
  method = c("grid", "random", "thinning"),
  cell_size = NULL,
  min_dist = NULL,
  seed = NULL
)
```

## Arguments

- coords:

  A data.frame with columns `x` and `y` containing site coordinates.

- n:

  Integer. Target number of sites to retain. If `NULL`, determined by
  `cell_size` or `min_dist`.

- method:

  Character. Subsampling method: `"grid"` (default), `"random"`, or
  `"thinning"`.

- cell_size:

  Numeric. Grid cell size for `method = "grid"`. One site retained per
  cell.

- min_dist:

  Numeric. Minimum distance between retained sites for
  `method = "thinning"`.

- seed:

  Integer. Random seed for reproducibility.

## Value

Integer vector of row indices to retain.

## Details

**Methods:**

- `"grid"`: Overlay a grid and retain one random site per cell.

- `"random"`: Simple random subsample of n sites.

- `"thinning"`: Iteratively remove sites until minimum distance is
  achieved.

## Examples

``` r
if (FALSE) { # \dontrun{
# Grid-based subsampling
keep <- subsample(coords, method = "grid", cell_size = 10)
sac <- spacc(species[keep, ], coords[keep, ])

# Minimum distance thinning
keep <- subsample(coords, method = "thinning", min_dist = 5)
} # }
```
