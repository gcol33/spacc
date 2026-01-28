# Phylogenetic Beta Diversity Accumulation

Compute spatial accumulation of phylogenetic beta diversity, partitioned
into turnover and nestedness components. Measures how evolutionary
composition changes as sites are accumulated spatially.

## Usage

``` r
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

- tree:

  A phylogenetic tree of class `phylo` (from ape), or a pairwise
  phylogenetic distance matrix.

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
`beta_type = "phylogenetic"`.

## Details

Phylogenetic beta diversity quantifies evolutionary turnover across
space. The PhyloSor index (phylogenetic Sorensen) is used: the fraction
of branch length shared between two communities relative to total branch
length. Partitioned into replacement (turnover) and loss (nestedness)
components.

## References

Baselga, A. (2010). Partitioning the turnover and nestedness components
of beta diversity. Global Ecology and Biogeography, 19, 134-143.

Chao, A., Chiu, C.H., Villeger, S., et al. (2023). Rarefaction and
extrapolation with beta diversity under a framework of Hill numbers: the
iNEXT.beta3D standardization. Ecological Monographs, 93, e1588.

## See also

[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md),
[`spaccBetaFunc()`](https://gillescolling.com/spacc/reference/spaccBetaFunc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ape)
tree <- rtree(20)
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 20, 1, 0.3), nrow = 50)
colnames(species) <- tree$tip.label

beta_phylo <- spaccBetaPhylo(species, coords, tree)
plot(beta_phylo)
} # }
```
