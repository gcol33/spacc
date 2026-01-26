# Parallel Radius-Order Accumulation

Parallel Radius-Order Accumulation

## Usage

``` r
cpp_radius_parallel(
  species_pa,
  dist_mat,
  n_seeds,
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

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

Integer matrix of accumulation curves
