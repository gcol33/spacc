# Parallel kNN Accumulation with Explicit Seeds

Run kNN accumulation from specified starting points in parallel.

## Usage

``` r
cpp_knn_parallel_seeds(
  species_pa,
  dist_mat,
  seeds,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Numeric matrix of pairwise distances

- seeds:

  Integer vector of starting point indices (0-based)

- n_cores:

  Number of cores to use

- progress:

  Show progress (currently ignored in C++)

## Value

Integer matrix (n_seeds x n_sites) of accumulation curves
