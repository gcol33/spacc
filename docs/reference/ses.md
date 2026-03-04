# Standardized Effect Size (SES) via Null Models

Compute standardized effect sizes by comparing observed diversity
metrics against a null distribution from randomized communities.
Supports multiple null model algorithms and works with most spacc output
classes.

## Usage

``` r
ses(
  x,
  species,
  coords = NULL,
  metric = NULL,
  null_model = c("frequency", "richness", "both", "curveball", "torus", "spatial_swap"),
  n_perm = 999L,
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  A spacc output object. Supported classes: `spacc`, `spacc_hill`,
  `spacc_phylo`, `spacc_func`, `spacc_beta`, `spacc_metrics`,
  `spacc_alpha`, `spacc_partition`.

- species:

  A site-by-species matrix (required). The species matrix used to
  produce `x`. Needed because spacc objects do not store the raw species
  matrix.

- coords:

  Optional data.frame with columns `x` and `y`. Required if the original
  analysis used coordinates. If `x` contains stored coordinates, they
  will be used automatically.

- metric:

  Character or `NULL`. For multi-metric objects (e.g., `spacc_hill` with
  multiple q), specify which metric to extract. If `NULL`, uses the
  first/default metric.

- null_model:

  Character. Null model algorithm:

  - `"frequency"`: Shuffle species columns independently (maintains
    column totals = species frequency)

  - `"richness"`: Shuffle species rows independently (maintains row
    totals = site richness)

  - `"both"`: Independent swap algorithm maintaining both row and column
    totals (Gotelli 2000)

  - `"curveball"`: Curveball algorithm for efficient swap (Strona et al.
    2014)

  - `"torus"`: Toroidal shift preserving spatial autocorrelation. Shifts
    all coordinates by a random offset and reassigns species to shifted
    sites. Requires `coords`.

  - `"spatial_swap"`: Independent swap restricted to spatially proximate
    site pairs. Preserves both marginals while respecting spatial
    structure. Requires `coords`.

- n_perm:

  Integer. Number of permutations. Default 999.

- parallel:

  Logical. Use parallel processing for the underlying analysis? Default
  `TRUE`.

- n_cores:

  Integer. Number of cores. Default `NULL` (auto-detect).

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed for reproducibility.

## Value

An object of class `spacc_ses` containing:

- observed:

  Numeric vector of observed metric values

- null_mean:

  Mean of null distribution

- null_sd:

  Standard deviation of null distribution

- ses:

  Standardized effect size: (observed - null_mean) / null_sd

- p_value:

  Two-tailed p-value

- n_perm:

  Number of permutations

- null_model:

  Null model algorithm used

- metric:

  Metric name

- input_class:

  Class of input object

## Details

SES is computed as: \$\$SES = \frac{observed -
\bar{null}}{sd\_{null}}\$\$

A two-tailed p-value is calculated as the proportion of null values at
least as extreme as the observed value: \$\$p = \frac{2 \cdot \min(r,
n\_{perm} + 1 - r)}{n\_{perm} + 1}\$\$ where \\r\\ is the rank of the
observed value among null values.

**Null model algorithms:**

- `"frequency"`: Tests whether species composition matters given
  observed species frequencies

- `"richness"`: Tests whether species identity matters given observed
  site richness

- `"both"`: Maintains both marginal totals; tests non-random species
  co-occurrence patterns

- `"curveball"`: Efficient alternative to `"both"` with proven uniform
  sampling properties

## References

Gotelli, N.J. (2000). Null model analysis of species co-occurrence
patterns. Ecology, 81, 2606-2621.

Strona, G., Nappo, D., Boccacci, F., Fattorini, S. & San-Miguel-Ayanz,
J. (2014). A fast and unbiased procedure to randomize ecological binary
matrices with fixed row and column totals. Nature Communications, 5,
4114.

## See also

[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md),
[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md),
[`spaccMetrics()`](https://gillescolling.com/spacc/reference/spaccMetrics.md)

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(20), y = runif(20))
species <- matrix(rbinom(20 * 15, 1, 0.3), nrow = 20)

sac <- spacc(species, coords, n_seeds = 10)
result <- ses(sac, species, n_perm = 19)
print(result)
# }
```
