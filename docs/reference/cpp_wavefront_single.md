# Single Wavefront Accumulation

Accumulate sites within expanding radius: r = r0 + step \* dr Models
invasion spread from introduction point.

## Usage

``` r
cpp_wavefront_single(species_pa, dist_mat, seed, r0, dr, n_steps)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Numeric matrix of pairwise distances

- seed:

  Starting site index (0-based)

- r0:

  Initial radius

- dr:

  Radius increment per step

- n_steps:

  Number of expansion steps

## Value

List with: species (counts), radius (values), n_sites (included)
