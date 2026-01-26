# Fast Distance Matrix

Compute pairwise distances between points using Euclidean or Haversine
formula.

## Usage

``` r
cpp_distance_matrix(x, y, method = "euclidean")
```

## Arguments

- x:

  Numeric vector of x coordinates (or longitude for haversine)

- y:

  Numeric vector of y coordinates (or latitude for haversine)

- method:

  Distance method: "euclidean" or "haversine"

## Value

Symmetric n x n distance matrix
