# Spatial Endemism Accumulation

Compute the number of endemic species (species found only within the
accumulated area) as sites are added spatially. Complements the standard
SAC by tracking species unique to each spatial extent.

## Usage

``` r
spaccEndemism(
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

  A site-by-species matrix (presence/absence or abundance).

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

An object of class `spacc_endemism` containing:

- richness:

  Matrix of cumulative richness (n_seeds x n_sites)

- endemism:

  Matrix of endemic species count (n_seeds x n_sites)

- coords, n_seeds, n_sites, method:

  Parameters used

## Details

At each accumulation step k, an endemic species is one that is present
in the accumulated sites (1..k) but absent from all remaining unvisited
sites (k+1..n). This tracks how many species are unique to the area
sampled so far.

The endemism curve typically starts low (few endemics at small areas),
increases as the region grows, and eventually equals total richness when
all sites are included.

## References

Kier, G., Kreft, H., Lee, T.M., et al. (2009). A global assessment of
endemism and species richness across island and mainland regions.
Proceedings of the National Academy of Sciences, 106, 9322-9327.

May, F., Gerstner, K., McGlinn, D.J., et al. (2018). mobsim: an R
package for the simulation and measurement of biodiversity across
spatial scales. Methods in Ecology and Evolution, 9, 1401-1408.

## See also

[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md),
[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md)

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

end <- spaccEndemism(species, coords, n_seeds = 30)
plot(end)
} # }
```
