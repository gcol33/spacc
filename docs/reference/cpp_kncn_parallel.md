# Parallel kNCN Accumulation

Run kNCN accumulation from multiple random starting points in parallel.

## Usage

``` r
cpp_kncn_parallel(species_pa, x, y, n_seeds, n_cores = 1L, progress = FALSE)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- x:

  Numeric vector of x coordinates

- y:

  Numeric vector of y coordinates

- n_seeds:

  Number of random starting points

- n_cores:

  Number of cores to use

- progress:

  Show progress (currently ignored)

## Value

Integer matrix (n_seeds x n_sites) of accumulation curves
