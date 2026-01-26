# Coverage-Based Spatial Rarefaction

Compute spatial accumulation curves with sample coverage tracking.
Allows standardization by completeness (coverage) rather than sample
size, following Chao & Jost (2012).

## Usage

``` r
spaccCoverage(
  x,
  coords,
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

  A site-by-species matrix with abundance data.

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method. Default `"knn"`.

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

## Value

An object of class `spacc_coverage` containing:

- richness:

  Matrix of species richness (n_seeds x n_sites)

- individuals:

  Matrix of individual counts

- coverage:

  Matrix of coverage estimates

- coords, n_seeds, n_sites, method:

  Parameters used

## Details

Sample coverage (Chao & Jost 2012) estimates the proportion of the total
community abundance represented by observed species. It provides a
measure of sampling completeness that is independent of sample size.

Coverage-based rarefaction allows fair comparison of diversity across
communities with different abundances by standardizing to equal
completeness rather than equal sample size.

## References

Chao, A. & Jost, L. (2012). Coverage-based rarefaction and
extrapolation: standardizing samples by completeness rather than size.
Ecology, 93, 2533-2547.

## See also

`iNEXT::iNEXT()` for coverage-based rarefaction without spatial
structure

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 30, 2), nrow = 50)

cov <- spaccCoverage(species, coords)
plot(cov)

# Interpolate richness at 90% and 95% coverage
interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
} # }
```
