# Spatial Hill Numbers at Standardized Coverage

Compute spatial accumulation of Hill numbers alongside sample coverage,
enabling standardized comparison at equal completeness levels
(iNEXT-style analysis applied spatially).

## Usage

``` r
spaccHillCoverage(
  x,
  coords,
  q = c(0, 1, 2),
  target_coverage = NULL,
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

  A site-by-species matrix (rows = sites, cols = species) with abundance
  data.

- coords:

  A data.frame with columns `x` and `y` containing site coordinates, or
  a `spacc_dist` object from
  [`distances()`](https://gillescolling.com/spacc/reference/distances.md).

- q:

  Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.

- target_coverage:

  Numeric vector. Target coverage levels for standardization. Default
  `NULL` (no standardization). Values in (0, 1).

- n_seeds:

  Integer. Number of random starting points. Default 50.

- distance:

  Character. Distance method: `"euclidean"` or `"haversine"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores. Default `NULL` uses all minus one.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed for reproducibility.

## Value

An object of class `spacc_hill_coverage` containing:

- hills:

  Named list of n_seeds x n_sites matrices (one per q)

- coverage:

  n_seeds x n_sites matrix of Chao-Jost coverage

- standardized:

  List of per-q values interpolated at target coverage, or `NULL` if
  `target_coverage` is `NULL`

- q:

  Vector of q values

- target_coverage:

  Target coverage levels used

- coords:

  Original coordinates

- n_seeds:

  Number of seeds

- n_sites:

  Number of sites

- n_species:

  Total species

## Details

Combines spatial kNN accumulation with simultaneous Hill number and
coverage computation. At each accumulation step, both the Chao-Jost
sample coverage and Hill numbers for all requested orders are
calculated.

When `target_coverage` is specified, Hill numbers are interpolated at
those coverage levels using the existing `interpolate_at_coverage()` C++
function, enabling fair comparison between sites with different sampling
completeness.

## References

Chao, A. & Jost, L. (2012). Coverage-based rarefaction and
extrapolation: standardizing samples by completeness rather than size.
Ecology, 93, 2533-2547.

Chao, A., Gotelli, N.J., Hsieh, T.C., Sander, E.L., Ma, K.H., Colwell,
R.K. & Ellison, A.M. (2014). Rarefaction and extrapolation with Hill
numbers. Ecological Monographs, 84, 45-67.

## See also

[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md)
for Hill accumulation without coverage,
[`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md)
for coverage accumulation without Hill numbers

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(40), y = runif(40))
species <- matrix(rpois(40 * 20, 2), nrow = 40)

hc <- spaccHillCoverage(species, coords, n_seeds = 10, progress = FALSE)
plot(hc)

# Standardize at 90% coverage
hc2 <- spaccHillCoverage(species, coords, target_coverage = 0.9,
                          n_seeds = 10, progress = FALSE)
summary(hc2)
# }
```
