# Spatial Accumulation with Hill Numbers

Compute spatial species accumulation curves using Hill numbers
(effective number of species) instead of raw richness. Hill numbers
unify diversity measures: q=0 is richness, q=1 is exponential Shannon,
q=2 is inverse Simpson.

## Usage

``` r
spaccHill(
  x,
  coords,
  q = c(0, 1, 2),
  n_seeds = 50L,
  method = "knn",
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (rows = sites, cols = species) with
  presence/absence (0/1) or abundance data.

- coords:

  A data.frame with columns `x` and `y` containing site coordinates, or
  a `spacc_dist` object from
  [`distances()`](https://gillescolling.com/spacc/reference/distances.md).

- q:

  Numeric vector. Orders of diversity to compute. Default `c(0, 1, 2)`.

  - q = 0: Species richness

  - q = 1: Exponential of Shannon entropy (effective common species)

  - q = 2: Inverse Simpson (effective dominant species)

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method: `"knn"` (default).

- distance:

  Character. Distance method: `"euclidean"` or `"haversine"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores. Default `NULL` uses `detectCores() - 1`.

- progress:

  Logical. Show progress bar? Default `TRUE`.

- seed:

  Integer. Random seed for reproducibility.

## Value

An object of class `spacc_hill` containing:

- curves:

  Named list of matrices, one per q value (n_seeds x n_sites)

- q:

  Vector of q values used

- coords:

  Original coordinates

- n_seeds:

  Number of seeds

- n_sites:

  Number of sites

- n_species:

  Total species

- method:

  Method used

## Details

Hill numbers (Chao et al. 2014) provide a unified framework for
diversity measurement. Unlike raw richness (q=0), higher-order Hill
numbers (q=1, q=2) down-weight rare species, providing different
perspectives on diversity.

The spatial accumulation of Hill numbers can reveal scale-dependent
diversity patterns missed by richness alone.

## References

Chao, A., Gotelli, N.J., Hsieh, T.C., Sander, E.L., Ma, K.H., Colwell,
R.K. & Ellison, A.M. (2014). Rarefaction and extrapolation with Hill
numbers: a framework for sampling and estimation in species diversity
studies. Ecological Monographs, 84, 45-67.

## See also

[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) for
richness-only accumulation, `iNEXT::iNEXT()` for non-spatial Hill number
rarefaction

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare diversity at different orders
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 30, 2), nrow = 50)

hill <- spaccHill(species, coords, q = c(0, 1, 2))
plot(hill)

# Extract summary at final site
summary(hill)
} # }
```
