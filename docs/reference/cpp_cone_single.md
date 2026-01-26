# Single Directional Cone Accumulation

Expand in a direction with angular constraint.

## Usage

``` r
cpp_cone_single(species_pa, x, y, seed, angle, width)
```

## Arguments

- species_pa:

  Integer matrix (sites x species)

- x:

  Numeric vector of x coordinates

- y:

  Numeric vector of y coordinates

- seed:

  Starting site index (0-based)

- angle:

  Direction angle in radians (0 = east, pi/2 = north)

- width:

  Cone half-width in radians

## Value

Integer vector of cumulative species counts
