# Parallel kNCN Metrics Accumulation

Run kNCN accumulation from each site as its own starting point.

## Usage

``` r
cpp_kncn_metrics_parallel(
  species_pa,
  x_coords,
  y_coords,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- x_coords:

  Numeric vector of x coordinates

- y_coords:

  Numeric vector of y coordinates

- n_cores:

  Number of cores to use

- progress:

  Show progress (currently ignored in C++)

## Value

Integer matrix (n_sites x n_sites) of accumulation curves
