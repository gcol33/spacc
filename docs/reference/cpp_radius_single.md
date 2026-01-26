# Single Radius-Order Accumulation

Accumulate sites in order of distance from seed (simpler version).

## Usage

``` r
cpp_radius_single(species_pa, dist_mat, seed)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Numeric matrix of pairwise distances

- seed:

  Starting site index (0-based)

## Value

Integer vector of cumulative species counts
