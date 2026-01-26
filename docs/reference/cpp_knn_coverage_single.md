# Single kNN Accumulation with Coverage Tracking

Single kNN Accumulation with Coverage Tracking

## Usage

``` r
cpp_knn_coverage_single(species_mat, dist_mat, seed)
```

## Arguments

- species_mat:

  Integer matrix (sites x species) with abundances

- dist_mat:

  Numeric matrix of pairwise distances

- seed:

  Starting site index (0-based)

## Value

List with richness, individuals, coverage vectors
