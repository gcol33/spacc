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
  method = c("knn", "nn_walk"),
  distance = c("euclidean", "haversine"),
  coverage = c("chao", "chiu"),
  n_units = NULL,
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  map = FALSE,
  focal_points = NULL,
  focal_domain = NULL
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

  Character. Accumulation method: canonical fixed-focus `"knn"` or
  `"nn_walk"`. Default `"knn"`.

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

- n_units:

  Total number of sampling units in the study region, when the sites are
  drawn without replacement from a finite set of units (atlas cells, a
  plot network). Requires `coverage = "chiu"`. `NULL` (default) or `Inf`
  assumes sampling with replacement. See Details.

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

- focal_points:

  Optional focal points for `method = "knn"`. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

- focal_domain:

  Optional polygonal focal domain. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

## Value

An object of class `spacc_coverage` containing:

- richness:

  Matrix of species richness (n_seeds x n_sites)

- individuals:

  Matrix of individual counts

- coverage:

  Matrix of coverage estimates

- pooled:

  Species abundances pooled over all sites

- coverage_type:

  Coverage estimator used (`"chao"` or `"chiu"`)

- n_units:

  Number of units in the finite population, or `NULL`

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
It writes the abundance of undetected species as the Chao2 estimate of
their number, \\\frac{t-1}{t} \frac{Q_1^2}{2 Q_2}\\ after \\t\\ sites,
times the mean abundance of a unique species, \\G_1 / Q_1\\.

With `n_units = T`, the \\t\\ sites are a fraction \\w = t / T\\ of a
finite set of units drawn without replacement, and the number of
undetected species is the MoRE term of `chao2(x, n_units = )` (Chiu
2023, Methods in Ecology and Evolution, Eq. 5), \$\$\frac{t-1}{t}
\frac{(1-w) Q_1^2}{Q_2 + w (Q_1 - Q_2)^+},\$\$ so richness and coverage
of the same sample rest on one estimate of the undetected species. It
vanishes at \\w = 1\\: the coverage reaches 1 once the accumulation has
visited every unit. On simulated grids its error is lower than that of
the with-replacement form from \\w = 0.5\\ on, and several-fold lower
near a census; it is about equal at \\w = 0.3\\ and larger at \\w =
0.1\\ on a sparse grid
([`vignette("extrapolation")`](https://gillescolling.com/spacc/articles/extrapolation.md)).

Coverage-based rarefaction allows fair comparison of diversity across
communities with different abundances by standardizing to equal
completeness rather than equal sample size.

## References

Chao, A. & Jost, L. (2012). Coverage-based rarefaction and
extrapolation: standardizing samples by completeness rather than size.
Ecology, 93, 2533-2547.

Chiu, C.-H. (2023). A species richness estimator for sample-based
incidence data sampled without replacement. Methods in Ecology and
Evolution, 14, 2482-2493.
[doi:10.1111/2041-210X.14146](https://doi.org/10.1111/2041-210X.14146)

Chiu, C.-H. (2023). A sample-based estimator for sample coverage.
Ecology, 104, e4099.

## See also

[`iNEXT::iNEXT()`](https://rdrr.io/pkg/iNEXT/man/iNEXT.html) for
coverage-based rarefaction without spatial structure

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 30, 2), nrow = 50)

cov <- spaccCoverage(species, coords)
plot(cov)

# Sample-based coverage (recommended for spatial data)
cov_chiu <- spaccCoverage(species, coords, coverage = "chiu")

# The 50 sites are 80% of a 62-cell grid
cov_grid <- spaccCoverage(species, coords, coverage = "chiu", n_units = 62)

# Interpolate richness at 90% and 95% coverage
interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
# }
```
