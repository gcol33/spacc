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
  rho = 0,
  n_seeds = 50L,
  method = c("knn", "nn_walk"),
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  map = FALSE,
  focal_points = NULL,
  focal_domain = NULL
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

  - `"evoheritage"`: EvoHeritage (Rosindell et al. 2024), PD with
    attrition of evolutionary heritage at rate `rho` (requires tree)

  - `"rao"`: Rao's quadratic entropy (abundance-weighted mean pairwise
    phylogenetic distance). Pass abundance data for weighting; with
    presence/absence it reduces to the equal-weight form.

- rho:

  Non-negative numeric vector. Attrition rate(s) for `"evoheritage"`,
  per unit of branch length. `rho = 0` gives Faith's PD. Several values
  give one curve each, named `evoheritage_rho<value>`.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method: `"knn"` or `"nn_walk"`.

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

- focal_points:

  Optional focal points for `method = "knn"`. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

- focal_domain:

  Optional polygonal focal domain. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

## Value

An object of class `spacc_phylo` containing:

- curves:

  Named list of matrices, one per metric (n_seeds x n_sites)

- metric:

  Metrics computed

- tree:

  The tree or distance matrix used

- coords, n_seeds, n_sites, method:

  Parameters used

## Details

Phylogenetic diversity metrics incorporate evolutionary relationships:

- **MPD (Mean Pairwise Distance)**: Average phylogenetic distance
  between all pairs of species. Sensitive to tree-wide patterns.

- **MNTD (Mean Nearest Taxon Distance)**: Average distance to closest
  relative. Sensitive to terminal clustering.

- **PD (Faith's Phylogenetic Diversity)**: Total branch length
  connecting species to the root. Requires full tree object.

- **EvoHeritage**: Features accumulate along each branch and are lost at
  rate `rho`, so a branch of length `L` passes on a fraction
  `exp(-rho * L)` of what entered it. A feature counts once if it
  survives to at least one present species. The value is in standardised
  units (the heritage at the end of a single branch of unit length), so
  `rho = 0` gives PD and large `rho` tends to species richness. Scanning
  `rho` places a community between its richness and PD curves.

Tree-based metrics match the columns of `x` to `tree$tip.label` by name;
the tree may contain species absent from `x`.

## References

Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity.
Biological Conservation, 61, 1-10.

Rosindell, J., Manson, K., Gumbs, R., Pearse, W.D. & Steel, M. (2024).
Phylogenetic biodiversity metrics should account for both accumulation
and attrition of evolutionary heritage. Systematic Biology, 73, 158-182.
[doi:10.1093/sysbio/syad072](https://doi.org/10.1093/sysbio/syad072)

Webb, C.O. (2000). Exploring the phylogenetic structure of ecological
communities: an example for rain forest trees. American Naturalist, 156,
145-155.

## See also

[`picante::mpd()`](https://rdrr.io/pkg/picante/man/mpd.html),
[`picante::mntd()`](https://rdrr.io/pkg/picante/man/mntd.html),
[`picante::pd()`](https://rdrr.io/pkg/picante/man/pd.html)

## Examples

``` r
# \donttest{
if (requireNamespace("ape", quietly = TRUE)) {
  # Create random tree
  tree <- ape::rtree(30)

  coords <- data.frame(x = runif(50), y = runif(50))
  species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
  colnames(species) <- tree$tip.label

  phylo <- spaccPhylo(species, coords, tree, metric = c("mpd", "mntd"))
  plot(phylo)
}
# }
```
