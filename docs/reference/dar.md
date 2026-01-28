# Diversity-Area Relationship (DAR)

Extend the classic species-area relationship (SAR) to a diversity-area
relationship using Hill numbers of any order q. Instead of plotting
species richness vs. sites, this plots effective diversity vs.
cumulative area.

## Usage

``` r
dar(
  x,
  coords,
  q = c(0, 1, 2),
  n_seeds = 50L,
  method = "knn",
  area_method = c("convex_hull", "voronoi", "count"),
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

- q:

  Numeric vector. Diversity orders. Default `c(0, 1, 2)`.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method. Default `"knn"`.

- area_method:

  Character. How to estimate cumulative area: `"voronoi"` (Voronoi
  tessellation, requires sf), `"convex_hull"` (convex hull of
  accumulated sites, requires sf), or `"count"` (use site count as
  proxy, no dependencies). Default `"convex_hull"`.

- distance:

  Character. Distance method. Default `"euclidean"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed.

## Value

An object of class `spacc_dar` containing:

- hill:

  A `spacc_hill` object with diversity curves

- area:

  Matrix of cumulative areas (n_seeds x n_sites)

- q:

  Diversity orders used

- area_method:

  Method used for area estimation

## Details

The DAR (Ma, 2018) generalizes the SAR by replacing species richness
(q=0) with Hill numbers of any order. This reveals how different facets
of diversity scale with area:

- q=0 (richness) recovers the classic SAR

- q=1 (Shannon) shows how common species diversity scales

- q=2 (Simpson) shows how dominant species diversity scales

## References

Ma, Z.S. (2018). DAR (diversity-area relationship): extending classic
SAR for biodiversity and biogeography analyses. Ecology and Evolution,
8, 10023-10038.

Arrhenius, O. (1921). Species and area. Journal of Ecology, 9, 95-99.

## See also

[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md),
[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 30, 2), nrow = 50)

dar <- dar(species, coords, q = c(0, 1, 2))
plot(dar)
} # }
```
