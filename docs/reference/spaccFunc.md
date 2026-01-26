# Spatial Functional Diversity Accumulation

Compute spatial accumulation of functional diversity metrics based on
traits.

## Usage

``` r
spaccFunc(
  x,
  coords,
  traits,
  metric = c("fdis", "fric"),
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

  A site-by-species matrix (abundance data recommended).

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- traits:

  A species-by-traits matrix. Row names should match species (columns of
  x).

- metric:

  Character vector. Metrics to compute:

  - `"fdis"`: Functional Dispersion (mean distance to centroid)

  - `"fric"`: Functional Richness (convex hull volume approximation)

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method. Default `"knn"`.

- distance:

  Character. Site distance method: `"euclidean"` or `"haversine"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed.

## Value

An object of class `spacc_func` containing:

- curves:

  Named list of matrices, one per metric (n_seeds x n_sites)

- metric:

  Metrics computed

- coords, n_seeds, n_sites, method:

  Parameters used

## Details

Functional diversity metrics quantify trait space occupation:

- **FDis (Functional Dispersion)**: Abundance-weighted mean distance
  from the community centroid in trait space. Captures functional
  divergence.

- **FRic (Functional Richness)**: Volume of trait space occupied (convex
  hull). Requires more species than traits to compute.

## References

Laliberté, E. & Legendre, P. (2010). A distance-based framework for
measuring functional diversity from multiple traits. Ecology, 91,
299-305.

## See also

`FD::dbFD()` for comprehensive functional diversity analysis

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 20, 2), nrow = 50)

# Trait matrix (species x traits)
traits <- matrix(rnorm(20 * 3), nrow = 20)
rownames(traits) <- paste0("sp", 1:20)
colnames(species) <- rownames(traits)

func <- spaccFunc(species, coords, traits, metric = c("fdis", "fric"))
plot(func)
} # }
```
