# Single kNCN Accumulation Curve

k-Nearest Centroid Neighbor: after each step, recalculate the centroid
of all visited sites, then find the unvisited site nearest to that
centroid.

## Usage

``` r
cpp_kncn_single(species_pa, x, y, seed)
```

## Arguments

- species_pa:

  Integer matrix (sites x species), presence/absence

- x:

  Numeric vector of x coordinates

- y:

  Numeric vector of y coordinates

- seed:

  Starting site index (0-based)

## Value

Integer vector of cumulative species counts
