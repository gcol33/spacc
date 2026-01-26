# Single kNN Phylogenetic Diversity Accumulation

Single kNN Phylogenetic Diversity Accumulation

## Usage

``` r
cpp_phylo_knn_single(species_pa, site_dist_mat, phylo_dist_mat, seed, metrics)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- site_dist_mat:

  Site distance matrix

- phylo_dist_mat:

  Phylogenetic/functional distance matrix (species x species)

- seed:

  Starting site

- metrics:

  Character vector of metrics: "mpd", "mntd"

## Value

List with accumulation curves for each metric
