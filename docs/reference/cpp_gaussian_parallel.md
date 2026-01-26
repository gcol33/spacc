# Parallel Gaussian-weighted Accumulation

Parallel Gaussian-weighted Accumulation

## Usage

``` r
cpp_gaussian_parallel(
  species_pa,
  dist_mat,
  n_seeds,
  sigma,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Numeric matrix of pairwise distances

- n_seeds:

  Number of random starting points

- sigma:

  Gaussian bandwidth parameter

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

Integer matrix of accumulation curves
