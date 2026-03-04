# Convert spacc_metrics to sf

Convert metrics to an sf object for spatial analysis and integration
with the areaOfEffect package.

## Usage

``` r
as_sf(x, crs = NULL)
```

## Arguments

- x:

  A `spacc_metrics` object.

- crs:

  Coordinate reference system. Default `NULL` (no CRS). Use EPSG codes
  like `4326` for WGS84 or `32631` for UTM zone 31N.

## Value

An sf object with POINT geometries and metric columns.

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
metrics <- spaccMetrics(species, coords)
if (requireNamespace("sf", quietly = TRUE)) {
  metrics_sf <- as_sf(metrics)
}
# }
```
