# Parallel kNN Accumulation

Run kNN accumulation from multiple random starting points in parallel.

## Usage

``` r
cpp_knn_parallel(species_pa, dist_mat, n_seeds, n_cores = 1L, progress = FALSE)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Numeric matrix of pairwise distances

- n_seeds:

  Number of random starting points

- n_cores:

  Number of cores to use

- progress:

  Show progress (currently ignored in C++)

## Value

Integer matrix (n_seeds x n_sites) of accumulation curves
