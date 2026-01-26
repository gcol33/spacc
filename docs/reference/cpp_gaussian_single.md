# Single Gaussian-weighted Accumulation

Probabilistically select next site weighted by distance (Gaussian
decay). Small sigma approximates kNN, large sigma approximates random.

## Usage

``` r
cpp_gaussian_single(species_pa, dist_mat, seed, sigma)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- dist_mat:

  Numeric matrix of pairwise distances

- seed:

  Starting site index (0-based)

- sigma:

  Gaussian bandwidth parameter

## Value

Integer vector of cumulative species counts
