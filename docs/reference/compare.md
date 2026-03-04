# Compare Two Accumulation Curves

Test whether two species accumulation curves differ significantly.

## Usage

``` r
compare(
  x,
  y,
  method = c("permutation", "bootstrap", "auc"),
  normalize = FALSE,
  n_perm = 999L,
  ...
)
```

## Arguments

- x:

  A `spacc` object.

- y:

  A `spacc` object.

- method:

  Character. Comparison method: `"permutation"` (default),
  `"bootstrap"`, or `"auc"` (area under curve difference).

- normalize:

  Logical. If `TRUE`, each seed's curve is divided by its final value
  before computing AUC, so that curves are compared on a \[0, 1\] scale
  (proportion of total species). This compares the **shape** of
  accumulation rather than absolute species counts. Default `FALSE`.

- n_perm:

  Integer. Number of permutations/bootstrap replicates. Default 999.

- ...:

  Additional arguments passed to comparison methods.

## Value

An object of class `spacc_comp` containing:

- x_name, y_name:

  Names of compared objects

- auc_diff:

  Difference in area under curve

- p_value:

  P-value from permutation test

- saturation_diff:

  Difference in saturation points

- method:

  Comparison method used

- normalized:

  Whether curves were normalized before comparison

## References

Colwell, R.K., Mao, C.X. & Chang, J. (2004). Interpolating,
extrapolating, and comparing incidence-based species accumulation
curves. Ecology, 85, 2717-2727.

Gotelli, N.J. & Colwell, R.K. (2011). Estimating species richness. In:
Biological Diversity: Frontiers in Measurement and Assessment, pp.
39-54. Oxford University Press.

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
sp_a <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
sp_b <- matrix(rbinom(50 * 30, 1, 0.5), nrow = 50)

sac_a <- spacc(sp_a, coords, n_seeds = 10)
sac_b <- spacc(sp_b, coords, n_seeds = 10)
comp <- compare(sac_a, sac_b)
print(comp)
# }
```
