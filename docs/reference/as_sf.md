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
if (FALSE) { # \dontrun{
metrics <- spaccMetrics(species, coords)

# Convert to sf with UTM projection
metrics_sf <- as_sf(metrics, crs = 32631)

# Use with areaOfEffect for spatial analysis
if (requireNamespace("areaOfEffect", quietly = TRUE)) {
  result <- areaOfEffect::aoe(metrics_sf, support = study_area)
}
} # }
```
