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
  seed = NULL,
  method = c("knn", "nn_walk"),
  focal_points = NULL,
  focal_domain = NULL,
  standardize = c("none", "coverage"),
  coverage = NULL
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

- method:

  Character. Accumulation method: `"knn"` or `"nn_walk"`.

- focal_points:

  Optional focal points for `method = "knn"`. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

- focal_domain:

  Optional polygonal focal domain. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

- standardize:

  Character. `"none"` (default) decomposes the observed samples;
  `"coverage"` standardises gamma and alpha to a common sample coverage
  at every step (iNEXT.beta3D, Chao et al. 2023). Needs integer
  abundances.

- coverage:

  Numeric in (0, 1\]. Fixed coverage for `standardize = "coverage"`, the
  same at every step. `NULL` (default) uses, at each step, the lower of
  the observed coverages of the pooled and joint samples.

## Value

An object of class `spacc_hill_beta` containing:

- gamma:

  Named list of n_seeds x n_sites matrices (one per q)

- alpha:

  Named list of n_seeds x n_sites matrices (one per q)

- beta:

  Named list of n_seeds x n_sites matrices (one per q)

- dissimilarity:

  With `standardize = "coverage"`: list of the four dissimilarities
  `1-CqN`, `1-UqN`, `1-VqN`, `1-SqN`, each a named list of n_seeds x
  n_sites matrices (one per q); `NULL` otherwise

- coverage:

  With `standardize = "coverage"`: n_seeds x n_sites matrix of the
  coverage used at each step

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

## Coverage standardisation

Observed beta depends on sampling completeness: gamma and alpha approach
their asymptotes at different rates, so the ratio moves along the curve
even where composition does not change. With `standardize = "coverage"`,
gamma is the Hill number of the pooled sample and alpha the Hill number
of the joint sample (every species-by-site cell treated as a unit)
divided by k (Chiu et al. 2014), both rarefied or extrapolated to the
same sample coverage C with the estimators of Chao et al. (2014); beta
is their ratio (Chao et al. 2023). By default C is, at each step, the
lower of the two observed coverages, so both values come from the
observed samples or their rarefaction. A fixed `coverage` makes the
curve comparable across steps and datasets, at the cost of extrapolation
where a step's sample is less complete; for q = 0, Chao et al. (2023)
advise against extrapolating beyond the coverage of twice the sample
size.

The standardised beta is transformed to four dissimilarities between 0
and 1 (Chao et al. 2023): `1-CqN` (Sorensen-type, the classic Sorensen
index at q = 0), `1-UqN` (Jaccard-type), `1-VqN` and `1-SqN`. They are
`NA` at the first step, where there is one site.

## References

Jost, L. (2007). Partitioning diversity into independent alpha and beta
components. Ecology, 88, 2427-2439.

Chiu, C.-H., Jost, L. & Chao, A. (2014). Phylogenetic beta diversity,
similarity, and differentiation measures based on Hill numbers.
Ecological Monographs, 84, 21-44.

Chao, A., Gotelli, N.J., Hsieh, T.C., et al. (2014). Rarefaction and
extrapolation with Hill numbers: a framework for sampling and estimation
in species diversity studies. Ecological Monographs, 84, 45-67.

Chao, A., Thorn, S., Chiu, C.-H., et al. (2023). Rarefaction and
extrapolation with beta diversity under a framework of Hill numbers: the
iNEXT.beta3D standardization. Ecological Monographs, 93, e1588.
[doi:10.1002/ecm.1588](https://doi.org/10.1002/ecm.1588)

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
