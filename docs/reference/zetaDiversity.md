# Zeta Diversity

Compute zeta diversity — the mean number of species shared across k
sites — for increasing orders of k. The zeta decline curve reveals
community assembly processes: exponential decline suggests stochastic
assembly, while power-law decline indicates niche-based assembly.

## Usage

``` r
zetaDiversity(
  x,
  coords,
  orders = 1:10,
  n_samples = 100L,
  method = c("knn", "random"),
  distance = c("euclidean", "haversine"),
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).
  Automatically binarized.

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- orders:

  Integer vector. Orders of zeta diversity to compute (number of sites
  in each combination). Default `1:10`.

- n_samples:

  Integer. Number of random combinations to sample per order. Default
  100.

- method:

  Character. Method for selecting k-site combinations: `"knn"`
  (spatially nearest sites) or `"random"` (random combinations). Default
  `"knn"`.

- distance:

  Character. Distance method: `"euclidean"` or `"haversine"`.

- seed:

  Integer. Random seed for reproducibility. Default `NULL`.

- progress:

  Logical. Show progress? Default `TRUE`.

## Value

An object of class `spacc_zeta` containing:

- zeta:

  Mean zeta values per order

- zeta_sd:

  Standard deviations per order

- orders:

  The k values

- n_samples:

  Number of samples per order

- ratio:

  Zeta ratio: zeta_k / zeta\_(k-1)

- decline:

  Data.frame with exponential and power-law fit statistics

- method:

  Method used

- n_sites:

  Number of sites

- n_species:

  Total species count

## Details

Zeta diversity of order k (\\\zeta_k\\) is the mean number of species
shared across k sites. Key properties:

- \\\zeta_1\\ = mean species richness per site

- \\\zeta_2\\ = mean number of species shared by any two sites

- \\\zeta_k\\ decreases monotonically with k

The zeta decline ratio (\\\zeta_k / \zeta\_{k-1}\\) is diagnostic:

- Constant ratio: exponential decline (stochastic assembly)

- Increasing ratio: power-law decline (deterministic/niche-based
  assembly)

The `knn` method selects spatially nearest k sites from each focal site,
which is ecologically meaningful for testing spatial turnover. The
`random` method samples random k-site combinations, providing a null
expectation.

## References

Hui, C. & McGeoch, M.A. (2014). Zeta diversity as a concept and metric
that unifies incidence-based biodiversity patterns. The American
Naturalist, 184, 684-694.

Latombe, G., McGeoch, M.A., Nipperess, D.A. & Hui, C. (2018). zetadiv:
an R package for computing compositional change across multiple sites,
assemblages or cases. bioRxiv, 324897.

## See also

[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
for pairwise beta diversity,
[`distanceDecay()`](https://gillescolling.com/spacc/reference/distanceDecay.md)
for distance-decay relationships

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(30), y = runif(30))
species <- matrix(rbinom(30 * 20, 1, 0.3), nrow = 30)
zeta <- zetaDiversity(species, coords, orders = 1:5, n_samples = 50)
plot(zeta)
# }
```
