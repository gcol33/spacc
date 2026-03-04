# Compute Distance Matrix

Pre-compute pairwise distances between sites for reuse across multiple
[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) calls.
Supports sf objects with accurate geodesic distances for global-scale
studies.

## Usage

``` r
distances(x, method = NULL, fun = NULL, which = NULL)
```

## Arguments

- x:

  Site locations. Can be:

  - A data.frame with columns `x` and `y`

  - An sf object with POINT geometries

  - An sfc_POINT object

- method:

  Character. Distance method:

  - `"euclidean"`: Euclidean distance (for projected coordinates)

  - `"haversine"`: Great-circle distance (for lat/lon, fast
    approximation

  - `"geodesic"`: Accurate ellipsoidal distance via sf/S2 (for global
    scale) Default is auto-detected from CRS when using sf objects.

- fun:

  Optional custom distance function. Must take two coordinate vectors
  (x, y) and return a distance matrix. Overrides `method`.

- which:

  For sf objects, column name containing the geometry. Default uses
  active geometry.

## Value

An object of class `spacc_dist` containing the distance matrix with
coordinates stored as an attribute.

## Details

For continental and global-scale studies, use sf objects with geographic
CRS (e.g., EPSG:4326). The function will automatically use accurate
geodesic distances via the S2 spherical geometry library.

For smaller study areas with projected coordinates (UTM, etc.),
Euclidean distance is appropriate and faster.

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
d <- distances(coords)

# Reuse for multiple analyses
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
sac <- spacc(species, d)
# }
```
