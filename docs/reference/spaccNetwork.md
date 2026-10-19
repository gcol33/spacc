# Spatial Accumulation of Interaction Diversity

Accumulate an interaction network over space: interaction richness and
interaction Hill numbers, sample coverage, the richness of the two
interacting levels (for example plants and pollinators) and connectance,
all along the same site orders.

## Usage

``` r
spaccNetwork(
  interactions,
  coords,
  q = c(0, 1, 2),
  n_seeds = 50L,
  method = c("knn", "nn_walk"),
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  focal_points = NULL,
  focal_domain = NULL
)
```

## Arguments

- interactions:

  A data frame of observed interactions with columns `site`, `lower`,
  `upper` and optionally `count` (interaction frequency, default 1 per
  row). `lower` and `upper` name the species of the two levels; an
  interaction is a `lower`-`upper` pair.

- coords:

  A data frame with columns `x` and `y`, one row per site. Sites in
  `interactions` are matched to `coords$site` if present, otherwise to
  the row names of `coords`, otherwise taken as row numbers.

- q:

  Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method: `"knn"` or `"nn_walk"`.

- distance:

  Character. `"euclidean"` or `"haversine"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed.

- focal_points:

  Optional focal points for `method = "knn"`. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

- focal_domain:

  Optional polygonal focal domain. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

## Value

An object of class `spacc_network` containing:

- interaction, lower, upper:

  Named lists (one per q) of n_seeds x n_sites matrices: Hill numbers of
  the interaction frequencies and of the summed frequencies of each
  level's species

- coverage:

  List with the n_seeds x n_sites sample coverage matrices of
  `interaction`, `lower` and `upper`

- connectance:

  n_seeds x n_sites matrix of interaction richness divided by the
  product of the two species richnesses

- q, coords, n_seeds, n_sites, n_interactions, n_lower, n_upper, method,
  orders:

  Parameters, dimensions and site orders

## Details

Each unique `lower`-`upper` pair is a unit of diversity and its
frequency its abundance, so the interaction Hill number of order q is
the ordinary Hill number of the interaction frequencies of the
accumulated sites, the effective number of equally frequent
interactions, and the sample coverage is the Chao & Jost (2012) coverage
of those frequencies (Chiu et al. 2023). The species curves of the two
levels use each species' summed interaction frequency. Connectance is
the realised share of the possible interactions between the species
detected so far.

## References

Chiu, C.-H., Chao, A., Vogel, S., Kriegel, P. & Thorn, S. (2023).
Quantifying and estimating ecological network diversity based on
incomplete sampling data. Philosophical Transactions of the Royal
Society B, 378, 20220183.
[doi:10.1098/rstb.2022.0183](https://doi.org/10.1098/rstb.2022.0183)

## See also

[`spaccHillCoverage()`](https://gillescolling.com/spacc/reference/spaccHillCoverage.md)
for the same quantities on species

## Examples

``` r
# \donttest{
set.seed(1)
coords <- data.frame(x = runif(30), y = runif(30))
interactions <- data.frame(
  site = sample(30, 400, replace = TRUE),
  lower = paste0("plant", sample(12, 400, replace = TRUE)),
  upper = paste0("bee", sample(20, 400, replace = TRUE))
)
net <- spaccNetwork(interactions, coords, n_seeds = 10)
net
plot(net)
# }
```
