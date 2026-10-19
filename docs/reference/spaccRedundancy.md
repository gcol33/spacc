# Spatial Phylogenetic or Functional Redundancy Accumulation

Track how much of the taxonomic diversity added along the spatial
accumulation is new phylogenetic or functional diversity, as the ratio
of the taxonomic Hill number to the phylogenetic or functional Hill
number of the accumulated community.

## Usage

``` r
spaccRedundancy(
  x,
  coords,
  tree = NULL,
  traits = NULL,
  q = c(0, 1, 2),
  n_seeds = 50L,
  method = c("knn", "nn_walk"),
  distance = c("euclidean", "haversine"),
  dist_method = c("euclidean", "gower"),
  normalize = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  focal_points = NULL,
  focal_domain = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (abundance or presence/absence).

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- tree:

  A phylogenetic tree of class `phylo` with branch lengths. Supply
  either `tree` or `traits`.

- traits:

  A species-by-traits matrix or data.frame, row names matching the
  columns of `x`.

- q:

  Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method: `"knn"` or `"nn_walk"`.

- distance:

  Character. Site distance method: `"euclidean"` or `"haversine"`.

- dist_method:

  Character. Trait distance for `traits`: `"euclidean"` on standardised
  traits (default) or `"gower"`.

- normalize:

  Logical. Rescale trait distances to a maximum of 1. Default `TRUE`.

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

An object of class `spacc_redundancy` containing:

- curves:

  Named list (one per q) of n_seeds x n_sites matrices of the ratio
  taxonomic / phylogenetic or functional Hill number

- taxonomic:

  The taxonomic Hill numbers, same layout

- dimension:

  The phylogenetic or functional Hill numbers, same layout

- type:

  `"phylogenetic"` or `"functional"`

- q, coords, n_seeds, n_sites, n_species, method, orders:

  Parameters and the site orders used

## Details

At each accumulation step the function computes, for every order q, the
taxonomic Hill number \\^qD\\ of the accumulated community and either
its phylogenetic Hill number \\^qD(T)\\ (Chao et al. 2010; effective
number of equally abundant, equally distinct lineages, with \\T\\ the
abundance-weighted mean root-to-tip depth) or its similarity-sensitive
functional Hill number \\^qD^Z\\ (Leinster & Cobbold 2012; with \\Z =
1 - d\\ from trait distances). The ratio \$\$R_q = {}^qD / {}^qD(T)
\quad \textrm{or} \quad {}^qD / {}^qD^Z\$\$ is at least 1, equal to 1
when species are maximally distinct (a star tree, or every pair of
species at the maximal trait distance), and grows as the species added
are close relatives or functionally similar to those already present. A
ratio rising along the curve says the area is adding species faster than
lineages or trait combinations.

## References

Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity measures
based on Hill numbers. Philosophical Transactions of the Royal Society
B, 365, 3599-3609.

Leinster, T. & Cobbold, C.A. (2012). Measuring diversity: the importance
of species similarity. Ecology, 93, 477-489.

Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of Hill
numbers to DNA-based diversity analyses. Molecular Ecology Resources,
19, 804-817.
[doi:10.1111/1755-0998.13014](https://doi.org/10.1111/1755-0998.13014)

## See also

[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md),
[`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md),
[`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md),
[`diversityProfile()`](https://gillescolling.com/spacc/reference/diversityProfile.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::rcoal(25)
  coords <- data.frame(x = runif(40), y = runif(40))
  species <- matrix(rpois(40 * 25, 1), nrow = 40,
                    dimnames = list(NULL, tree$tip.label))
  red <- spaccRedundancy(species, coords, tree = tree, n_seeds = 10)
  plot(red)
}
# }
```
