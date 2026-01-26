# Parallel Directional Cone Accumulation

Parallel Directional Cone Accumulation

## Usage

``` r
cpp_cone_parallel(
  species_pa,
  x,
  y,
  n_seeds,
  width = 0.785398,
  n_cores = 1L,
  progress = FALSE
)
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

- width:

  Cone half-width in radians (default pi/4)

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

Integer matrix of accumulation curves
