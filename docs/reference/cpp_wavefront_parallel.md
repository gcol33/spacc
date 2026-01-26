# Parallel Wavefront Accumulation

Parallel Wavefront Accumulation

## Usage

``` r
cpp_wavefront_parallel(
  species_pa,
  dist_mat,
  n_seeds,
  r0,
  dr,
  n_steps,
  n_cores = 1L,
  progress = FALSE
)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Distance matrix

- n_seeds:

  Number of random starting points

- r0:

  Initial radius

- dr:

  Radius increment per step

- n_steps:

  Number of expansion steps

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

List with curves matrix and radius vector
