# Distance-Decay Accumulation

Returns species count at each distance threshold from seed.

## Usage

``` r
cpp_distance_decay_single(species_pa, dist_mat, seed, breaks)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Distance matrix

- seed:

  Starting site

- breaks:

  Distance thresholds

## Value

Integer vector of cumulative species at each threshold
