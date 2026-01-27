# Spatial Phylogenetic Diversity Accumulation

Compute spatial accumulation of phylogenetic diversity metrics (MPD,
MNTD, PD).

## Usage

``` r
spaccPhylo(
  x,
  coords,
  tree,
  metric = c("mpd", "mntd"),
  n_seeds = 50L,
  method = "knn",
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  map = FALSE
)
```

## Arguments

- x:

  A site-by-species matrix.

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- tree:

  A phylogenetic tree of class `phylo` (from ape package), or a pairwise
  phylogenetic distance matrix.

- metric:

  Character vector. Metrics to compute:

  - `"mpd"`: Mean Pairwise Distance

  - `"mntd"`: Mean Nearest Taxon Distance

  - `"pd"`: Faith's Phylogenetic Diversity (requires tree, not distance
    matrix)

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method. Default `"knn"`.

- distance:

  Character. Site distance method: `"euclidean"` or `"haversine"`.

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
  per-site final values for spatial mapping. Enables
  [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) and
  `plot(type = "map")`. Default `FALSE`.

## Value

An object of class `spacc_phylo` containing:

- curves:

  Named list of matrices, one per metric (n_seeds x n_sites)

- metric:

  Metrics computed

- coords, n_seeds, n_sites, method:

  Parameters used

## Details

Phylogenetic diversity metrics incorporate evolutionary relationships:

- **MPD (Mean Pairwise Distance)**: Average phylogenetic distance
  between all pairs of species. Sensitive to tree-wide patterns.

- **MNTD (Mean Nearest Taxon Distance)**: Average distance to closest
  relative. Sensitive to terminal clustering.

- **PD (Faith's Phylogenetic Diversity)**: Total branch length
  connecting species. Requires full tree object.

## References

Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity.
Biological Conservation, 61, 1-10.

Webb, C.O. (2000). Exploring the phylogenetic structure of ecological
communities: an example for rain forest trees. American Naturalist, 156,
145-155.

## See also

[`picante::mpd()`](https://rdrr.io/pkg/picante/man/mpd.html),
[`picante::mntd()`](https://rdrr.io/pkg/picante/man/mntd.html),
[`picante::pd()`](https://rdrr.io/pkg/picante/man/pd.html)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ape)

# Create random tree
tree <- rtree(30)

coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
colnames(species) <- tree$tip.label

phylo <- spaccPhylo(species, coords, tree, metric = c("mpd", "mntd"))
plot(phylo)
} # }
```
