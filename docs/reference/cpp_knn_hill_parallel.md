# Parallel kNN Accumulation with Hill Numbers

Parallel kNN Accumulation with Hill Numbers

## Usage

``` r
cpp_knn_hill_parallel(
  species_mat,
  dist_mat,
  n_seeds,
  q_values,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_mat:

  Integer matrix (sites x species) with abundances

- dist_mat:

  Numeric matrix of pairwise distances

- n_seeds:

  Number of random starting points

- q_values:

  Vector of Hill number orders

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

List with curves for each q value
