# Spatial Hill Number Beta Diversity

Compute multisite beta diversity as gamma/alpha decomposition of Hill
numbers along the spatial accumulation curve (Jost 2007 framework).

## Usage

``` r
spaccHillBeta(
  x,
  coords,
  q = c(0, 1, 2),
  n_seeds = 50L,
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (rows = sites, cols = species).

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- q:

  Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- distance:

  Character. `"euclidean"` or `"haversine"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores. Default `NULL`.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed.

## Value

An object of class `spacc_hill_beta` containing:

- gamma:

  Named list of n_seeds x n_sites matrices (one per q)

- alpha:

  Named list of n_seeds x n_sites matrices (one per q)

- beta:

  Named list of n_seeds x n_sites matrices (one per q)

- q:

  Vector of q values

- coords:

  Original coordinates

- n_seeds, n_sites, n_species:

  Dimensions

## Details

At each accumulation step k, the function computes:

- **Gamma**: Hill number of the pooled community (all k sites combined)

- **Alpha**: Generalized mean of per-site Hill numbers (Jost's power
  mean)

- **Beta**: gamma / alpha (effective number of distinct communities)

Beta = 1 means all sites are identical; beta = k means all sites are
completely different. This provides the Hill-number analogue of the
Baselga-based
[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md).

## References

Jost, L. (2007). Partitioning diversity into independent alpha and beta
components. Ecology, 88, 2427-2439.

## See also

[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
for P/A-based Baselga partitioning,
[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md)
for Hill accumulation without beta decomposition

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(40), y = runif(40))
species <- matrix(rpois(40 * 20, 2), nrow = 40)

hb <- spaccHillBeta(species, coords, n_seeds = 10, progress = FALSE)
plot(hb)
# }
```
