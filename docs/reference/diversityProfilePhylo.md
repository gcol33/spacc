# Phylogenetic Diversity Profile

Compute phylogenetic Hill numbers (Chao et al. 2010) across a continuous
range of diversity orders (q), producing a phylogenetic diversity
profile.

## Usage

``` r
diversityProfilePhylo(
  x,
  tree,
  q = seq(0, 3, by = 0.1),
  type = c("both", "per_site", "regional"),
  coords = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (abundance data). Column names must match tip
  labels in the phylogeny.

- tree:

  An [`ape::phylo`](https://rdrr.io/pkg/ape/man/read.tree.html) object.
  Tips must include all species in `x`.

- q:

  Numeric vector. Orders of diversity. Default `seq(0, 3, by = 0.1)`.

- type:

  Character. What to compute: `"per_site"`, `"regional"`, or `"both"`
  (default).

- coords:

  Optional data.frame with `x` and `y` for spatial mapping.

## Value

An object of class `spacc_profile` with
`$profile_type = "phylogenetic"`.

## Details

Phylogenetic Hill numbers (Chao et al. 2010) weight branches by their
evolutionary distance. At q=0 this approximates normalized Faith's PD.
Higher q values increasingly emphasize common lineages.

## References

Chao, A., Chiu, C.H. & Jost, L. (2010). Phylogenetic diversity measures
based on Hill numbers. Philosophical Transactions of the Royal Society
B, 365, 3599-3609.

## See also

[`diversityProfile()`](https://gillescolling.com/spacc/reference/diversityProfile.md)
for taxonomic profiles,
[`diversityProfileFunc()`](https://gillescolling.com/spacc/reference/diversityProfileFunc.md)
for functional profiles

## Examples

``` r
# \donttest{
if (requireNamespace("ape", quietly = TRUE)) {
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  tree <- ape::rcoal(10, tip.label = paste0("sp", 1:10))
  prof <- diversityProfilePhylo(species, tree)
  print(prof)
}
# }
```
