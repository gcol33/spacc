# Spatial Beta Diversity Accumulation

Analyze how beta diversity changes as sites are accumulated spatially.
Partitions beta diversity into turnover (species replacement) and
nestedness (species loss) components following Baselga (2010).

## Usage

``` r
spaccBeta(
  x,
  coords,
  traits = NULL,
  tree = NULL,
  n_seeds = 50L,
  method = "knn",
  index = c("sorensen", "jaccard"),
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  map = FALSE
)

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

spaccBetaPhylo(
  x,
  coords,
  tree,
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

  Optional species-by-traits matrix or data.frame (row names matching
  species). When supplied, functional beta diversity is computed.

- tree:

  Optional phylogenetic tree of class `phylo`, or a pairwise
  phylogenetic distance matrix. When supplied, phylogenetic beta
  diversity is computed. Supply at most one of `traits` or `tree`.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method. Default `"knn"`.

- index:

  Character. Dissimilarity index: `"sorensen"` (default) or `"jaccard"`.

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

- map:

  Logical. If `TRUE`, run accumulation from every site as seed and store
  per-site final beta values for spatial mapping. Enables
  [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) and
  `plot(type = "map")`. Default `FALSE`.

## Value

An object of class `spacc_beta` containing:

- beta_total:

  Matrix of total beta diversity (n_seeds x n_sites-1)

- beta_turnover:

  Matrix of turnover component

- beta_nestedness:

  Matrix of nestedness component

- distance:

  Matrix of cumulative distances

- n_seeds, n_sites, method, index:

  Parameters used

## Details

At each step of spatial accumulation, beta diversity is calculated
between the accumulated species pool and the newly added site. This
reveals how species composition changes as you expand spatially.

**Interpretation:**

- High turnover: New sites contribute different species (replacement)

- High nestedness: New sites contribute subsets of existing species
  (loss)

The sum of turnover and nestedness equals total beta diversity.

Supplying `traits` computes functional beta diversity, weighting the
partition by the trait distinctiveness of the species exchanged between
the accumulated pool and each new site. Supplying `tree` computes
phylogenetic beta diversity (PhyloSor), weighting by shared branch
length. The taxonomic, functional, and phylogenetic variants all return
a `spacc_beta` object whose `beta_type` field records which was
computed. `map = TRUE` is supported for the taxonomic variant only.

## References

Baselga, A. (2010). Partitioning the turnover and nestedness components
of beta diversity. Global Ecology and Biogeography, 19, 134-143.

## See also

[`betapart::beta.pair()`](https://rdrr.io/pkg/betapart/man/beta.pair.html)
for pairwise beta diversity

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

beta <- spaccBeta(species, coords, n_seeds = 30)
plot(beta)

# Compare Sorensen vs Jaccard
beta_jac <- spaccBeta(species, coords, index = "jaccard")

# Functional beta diversity (pass a species-by-traits matrix)
traits <- matrix(rnorm(30 * 3), nrow = 30)
rownames(traits) <- colnames(species) <- paste0("sp", 1:30)
beta_func <- spaccBeta(species, coords, traits = traits)
# }
```
