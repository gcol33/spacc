# Single kNN Accumulation Curve

Traverse sites in nearest-neighbor order from a starting seed,
accumulating species counts.

## Usage

``` r
cpp_knn_single(species_pa, dist_mat, seed)
```

## Arguments

- species_pa:

  Integer matrix (sites x species), presence/absence

- dist_mat:

  Numeric matrix of pairwise distances

- seed:

  Starting site index (0-based)

## Value

Integer vector of cumulative species counts
