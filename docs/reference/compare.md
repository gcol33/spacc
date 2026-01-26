# Compare Two Accumulation Curves

Test whether two species accumulation curves differ significantly.

## Usage

``` r
compare(
  x,
  y,
  method = c("permutation", "bootstrap", "auc"),
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

## Examples

``` r
if (FALSE) { # \dontrun{
sac_native <- spacc(native_species, coords)
sac_alien <- spacc(alien_species, coords)

comp <- compare(sac_native, sac_alien)
print(comp)
plot(comp)
} # }
```
