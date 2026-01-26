# Parallel kNN Beta Diversity Accumulation

Parallel kNN Beta Diversity Accumulation

## Usage

``` r
cpp_beta_knn_parallel(
  species_pa,
  dist_mat,
  n_seeds,
  use_jaccard = FALSE,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Distance matrix

- n_seeds:

  Number of starting points

- use_jaccard:

  Use Jaccard instead of Sorensen

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

List with matrices for each beta component
