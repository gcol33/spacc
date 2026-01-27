# Spatial Beta Diversity Accumulation

Analyze how beta diversity changes as sites are accumulated spatially.
Partitions beta diversity into turnover (species replacement) and
nestedness (species loss) components following Baselga (2010).

## Usage

``` r
spaccBeta(
  x,
  coords,
  n_seeds = 50L,
  method = "knn",
  index = c("sorensen", "jaccard"),
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  map = FALSE
)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method. Default `"knn"`.

- index:

  Character. Dissimilarity index: `"sorensen"` (default) or `"jaccard"`.

- distance:

  Character. Distance method: `"euclidean"` or `"haversine"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed.

- map:

  Logical. If `TRUE`, run accumulation from every site as seed and store
  per-site final beta values for spatial mapping. Enables
  [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) and
  `plot(type = "map")`. Default `FALSE`.

## Value

An object of class `spacc_beta` containing:

- beta_total:

  Matrix of total beta diversity (n_seeds x n_sites-1)

- beta_turnover:

  Matrix of turnover component

- beta_nestedness:

  Matrix of nestedness component

- distance:

  Matrix of cumulative distances

- n_seeds, n_sites, method, index:

  Parameters used

## Details

At each step of spatial accumulation, beta diversity is calculated
between the accumulated species pool and the newly added site. This
reveals how species composition changes as you expand spatially.

**Interpretation:**

- High turnover: New sites contribute different species (replacement)

- High nestedness: New sites contribute subsets of existing species
  (loss)

The sum of turnover and nestedness equals total beta diversity.

## References

Baselga, A. (2010). Partitioning the turnover and nestedness components
of beta diversity. Global Ecology and Biogeography, 19, 134-143.

## See also

[`betapart::beta.pair()`](https://rdrr.io/pkg/betapart/man/beta.pair.html)
for pairwise beta diversity

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

beta <- spaccBeta(species, coords, n_seeds = 30)
plot(beta)

# Compare Sorensen vs Jaccard
beta_jac <- spaccBeta(species, coords, index = "jaccard")
} # }
```
