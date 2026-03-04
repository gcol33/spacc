# Spatial Variance Partitioning with MEMs

Partition diversity variation into spatial and non-spatial components
using forward selection of Moran's Eigenvector Maps.

## Usage

``` r
spatialPartition(x, mem, metric = NULL, forward = TRUE, alpha = 0.05)
```

## Arguments

- x:

  A spacc result object (e.g., `spacc_metrics`, `spacc_alpha`,
  `spacc_hill`) or a numeric vector of per-site diversity values.

- mem:

  A `spacc_mem` object from
  [`spatialEigenvectors()`](https://gillescolling.com/spacc/reference/spatialEigenvectors.md).

- metric:

  Character. Which metric to extract from `x` (passed to internal
  extraction). Default `NULL` uses the first available.

- forward:

  Logical. Use forward selection? Default `TRUE`. If `FALSE`, all MEMs
  are included.

- alpha:

  Numeric. Significance threshold for forward selection. Default 0.05.

## Value

An object of class `spacc_mem_partition` containing:

- r_squared_spatial:

  R-squared of the spatial model

- r_squared_total:

  Total R-squared (same as spatial here)

- selected_mems:

  Names of selected MEM vectors

- n_selected:

  Number of selected MEMs

- anova_table:

  ANOVA table from the final model

- coefficients:

  Model coefficients

## Details

Forward selection of MEMs proceeds by adding the MEM that most improves
the model AIC at each step, stopping when no MEM improves AIC by more
than 2 units or when p \> alpha.

## See also

[`spatialEigenvectors()`](https://gillescolling.com/spacc/reference/spatialEigenvectors.md)
for computing MEMs

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(30), y = runif(30))
species <- matrix(rpois(30 * 15, 2), nrow = 30)

mem <- spatialEigenvectors(coords)
alpha <- alphaDiversity(species, q = 0)
part <- spatialPartition(alpha$q0, mem)
print(part)
# }
```
