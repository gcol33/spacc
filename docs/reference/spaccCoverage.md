# Coverage-Based Spatial Rarefaction

Compute spatial accumulation curves with sample coverage tracking.
Allows standardization by completeness (coverage) rather than sample
size, following Chao & Jost (2012) or the sample-based estimator of Chiu
(2023).

## Usage

``` r
spaccCoverage(
  x,
  coords,
  n_seeds = 50L,
  method = "knn",
  distance = c("euclidean", "haversine"),
  coverage = c("chao", "chiu"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  map = FALSE
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

- coverage:

  Character. Coverage estimator to use: `"chao"` (default) for the
  individual-based Chao & Jost (2012) estimator using
  singletons/doubletons, or `"chiu"` for the sample-based Chiu (2023)
  estimator using incidence frequency counts (Q1/Q2) and unique-species
  abundance (G1). The Chiu estimator is recommended for spatially
  aggregated data where sampling units are plots/sites rather than
  independent individuals.

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
  per-site final coverage and richness for spatial mapping. Enables
  [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) and
  `plot(type = "map")`. Default `FALSE`.

## Value

An object of class `spacc_coverage` containing:

- richness:

  Matrix of species richness (n_seeds x n_sites)

- individuals:

  Matrix of individual counts

- coverage:

  Matrix of coverage estimates

- coverage_type:

  Coverage estimator used (`"chao"` or `"chiu"`)

- coords, n_seeds, n_sites, method:

  Parameters used

## Details

Sample coverage estimates the proportion of the total community
abundance represented by observed species. It provides a measure of
sampling completeness that is independent of sample size.

The **Chao-Jost (2012)** estimator counts singletons (f1) and doubletons
(f2) in the cumulative abundance vector. It assumes individuals are
sampled independently, which may not hold for plot-based spatial data.

The **Chiu (2023)** estimator uses incidence frequency counts instead:
Q1 (species in exactly 1 site), Q2 (species in exactly 2 sites), and G1
(total abundance of Q1 species). It gives near-unbiased coverage
estimates when organisms are spatially aggregated across sampling units.

Coverage-based rarefaction allows fair comparison of diversity across
communities with different abundances by standardizing to equal
completeness rather than equal sample size.

## References

Chao, A. & Jost, L. (2012). Coverage-based rarefaction and
extrapolation: standardizing samples by completeness rather than size.
Ecology, 93, 2533-2547.

Chiu, C.-H. (2023). A sample-based estimator for sample coverage.
Ecology, 104, e4099.

## See also

`iNEXT::iNEXT()` for coverage-based rarefaction without spatial
structure

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 30, 2), nrow = 50)

cov <- spaccCoverage(species, coords)
plot(cov)

# Sample-based coverage (recommended for spatial data)
cov_chiu <- spaccCoverage(species, coords, coverage = "chiu")

# Interpolate richness at 90% and 95% coverage
interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
# }
```
