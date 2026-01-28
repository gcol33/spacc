# Functional Beta Diversity Accumulation

Compute spatial accumulation of functional beta diversity, partitioned
into turnover and nestedness components. Measures how functional trait
space composition changes as sites are accumulated spatially.

## Usage

``` r
spaccBetaFunc(
  x,
  coords,
  traits,
  n_seeds = 50L,
  method = "knn",
  index = c("sorensen", "jaccard"),
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

- traits:

  A species-by-traits matrix. Row names should match species.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method. Default `"knn"`.

- index:

  Character. Dissimilarity index: `"sorensen"` (default) or `"jaccard"`.

- distance:

  Character. Distance method. Default `"euclidean"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed.

## Value

An object of class `spacc_beta` with additional attribute
`beta_type = "functional"`.

## Details

Functional beta diversity quantifies the turnover of functional traits
across space. At each accumulation step, beta is computed based on the
overlap of trait ranges (functional space) between the accumulated pool
and the newly added site.

## References

Baselga, A. (2012). The relationship between species replacement,
dissimilarity derived from nestedness, and nestedness. Global Ecology
and Biogeography, 21, 1223-1232.

Cardoso, P., Rigal, F. & Carvalho, J.C. (2015). BAT – Biodiversity
Assessment Tools. Methods in Ecology and Evolution, 6, 232-236.

## See also

[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md),
[`spaccBetaPhylo()`](https://gillescolling.com/spacc/reference/spaccBetaPhylo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 20, 1, 0.3), nrow = 50)
traits <- matrix(rnorm(20 * 3), nrow = 20)
rownames(traits) <- colnames(species) <- paste0("sp", 1:20)

beta_func <- spaccBetaFunc(species, coords, traits)
plot(beta_func)
} # }
```
