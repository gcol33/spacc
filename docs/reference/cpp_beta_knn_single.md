# Single kNN Beta Diversity Accumulation

Track how beta diversity (total, turnover, nestedness) changes as sites
are accumulated spatially using kNN method.

## Usage

``` r
cpp_beta_knn_single(species_pa, dist_mat, seed, use_jaccard = FALSE)
```

## Arguments

- species_pa:

  Integer matrix (sites x species) presence/absence

- dist_mat:

  Numeric matrix of pairwise distances

- seed:

  Starting site index (0-based)

- use_jaccard:

  If true, use Jaccard; otherwise Sorensen

## Value

List with beta_total, beta_turnover, beta_nestedness, distance vectors
