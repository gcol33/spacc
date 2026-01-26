# Parallel Functional Diversity Accumulation

Parallel Functional Diversity Accumulation

## Usage

``` r
cpp_func_knn_parallel(
  species_mat,
  site_dist_mat,
  traits,
  n_seeds,
  metrics,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_mat:

  Integer matrix (sites x species) abundances

- site_dist_mat:

  Site distance matrix

- traits:

  Trait matrix (species x traits)

- n_seeds:

  Number of random starting points

- metrics:

  Character vector of metrics to compute

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

List with matrices for each metric
