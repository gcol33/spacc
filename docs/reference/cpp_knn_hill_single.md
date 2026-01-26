# Single kNN Accumulation with Hill Numbers

Single kNN Accumulation with Hill Numbers

## Usage

``` r
cpp_knn_hill_single(species_mat, dist_mat, seed, q_values)
```

## Arguments

- species_mat:

  Integer matrix (sites x species) with abundances

- dist_mat:

  Numeric matrix of pairwise distances

- seed:

  Starting site index (0-based)

- q_values:

  Vector of Hill number orders to compute

## Value

Numeric matrix (length(q) x n_sites) of Hill numbers
