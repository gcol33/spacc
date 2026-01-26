# Parallel Phylogenetic Diversity Accumulation

Parallel Phylogenetic Diversity Accumulation

## Usage

``` r
cpp_phylo_knn_parallel(
  species_pa,
  site_dist_mat,
  phylo_dist_mat,
  n_seeds,
  metrics,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- site_dist_mat:

  Site distance matrix

- phylo_dist_mat:

  Phylogenetic distance matrix

- n_seeds:

  Number of starting points

- metrics:

  Metrics to compute

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

List with matrices for each metric
