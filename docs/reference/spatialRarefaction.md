# Spatially-Constrained Rarefaction

Compute expected species richness accounting for spatial autocorrelation
(Chiarucci et al. 2009). Uses distance-weighted sampling probabilities.

## Usage

``` r
spatialRarefaction(x, coords, n_perm = 100, bandwidth = NULL)
```

## Arguments

- x:

  A site-by-species matrix.

- coords:

  A data.frame with x and y columns.

- n_perm:

  Number of permutations. Default 100.

- bandwidth:

  Distance bandwidth for spatial weighting.

## Value

A data.frame with columns: sites, mean, sd, lower, upper

## References

Chiarucci, A., Bacaro, G., Rocchini, D. & Fattorini, L. (2009).
Discovering and rediscovering the sample-based rarefaction formula in
the ecological literature. Community Ecology, 10, 195-199.
