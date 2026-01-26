# Parallel kNN Metrics Accumulation

Run kNN accumulation from each site as its own starting point. Returns
one curve per site for extracting per-site metrics.

## Usage

``` r
cpp_knn_metrics_parallel(species_pa, dist_mat, n_cores = 1L, progress = FALSE)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Numeric matrix of pairwise distances

- n_cores:

  Number of cores to use

- progress:

  Show progress (currently ignored in C++)

## Value

Integer matrix (n_sites x n_sites) of accumulation curves
