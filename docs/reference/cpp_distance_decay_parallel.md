# Parallel Distance-Decay Accumulation

Parallel Distance-Decay Accumulation

## Usage

``` r
cpp_distance_decay_parallel(
  species_pa,
  dist_mat,
  n_seeds,
  breaks,
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

- breaks:

  Distance thresholds

- n_cores:

  Number of cores

- progress:

  Show progress

## Value

Integer matrix of species counts
