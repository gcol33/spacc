# Parallel Random Accumulation

Run random-order accumulation multiple times (null model).

## Usage

``` r
cpp_random_parallel(species_pa, n_seeds, n_cores = 1L, progress = FALSE)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- n_seeds:

  Number of random permutations

- n_cores:

  Number of cores to use

- progress:

  Show progress (currently ignored)

## Value

Integer matrix (n_seeds x n_sites) of accumulation curves
