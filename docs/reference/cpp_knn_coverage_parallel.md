# Parallel kNN Accumulation with Coverage

Parallel kNN Accumulation with Coverage

## Usage

``` r
cpp_knn_coverage_parallel(
  species_mat,
  dist_mat,
  n_seeds,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_mat:

  Integer matrix (sites x species)

- dist_mat:

  Distance matrix

- n_seeds:

  Number of starting points

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

List with richness, individuals, coverage matrices
