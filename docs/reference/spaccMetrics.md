# Per-Site Accumulation Metrics

Compute spatial accumulation metrics for each site as a starting point.
Useful for identifying sites with high or low accumulation rates,
visualizing spatial patterns in diversity, and understanding edge
effects.

## Usage

``` r
spaccMetrics(
  x,
  coords,
  metrics = c("slope_10", "half_richness", "auc"),
  method = c("knn", "kncn", "random"),
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE
)
```

## Arguments

- x:

  A site-by-species matrix (rows = sites, cols = species).

- coords:

  A data.frame with columns `x` and `y` containing site coordinates, or
  a `spacc_dist` object from
  [`distances()`](https://gillescolling.com/spacc/reference/distances.md).

- metrics:

  Character vector. Metrics to compute. Options include: `"slope_10"`
  (initial slope, first 10% of sites), `"slope_25"` (initial slope,
  first 25% of sites), `"half_richness"` (sites to reach 50% of total
  species), `"richness_50pct"` (alias for half_richness),
  `"richness_75pct"` (sites to reach 75% of species), `"richness_90pct"`
  (sites to reach 90% of species), `"auc"` (area under accumulation
  curve), `"final_richness"` (total species starting from this site).

- method:

  Character. Accumulation method: `"knn"`, `"kncn"`, `"random"`. Default
  `"knn"`.

- distance:

  Character. Distance method: `"euclidean"` or `"haversine"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores for parallel processing.

- progress:

  Logical. Show progress bar? Default `TRUE`.

## Value

An object of class `spacc_metrics` containing:

- metrics:

  Data frame with one row per site and columns for each metric

- coords:

  Original coordinates

- metric_names:

  Names of computed metrics

- n_sites:

  Number of sites

- n_species:

  Total species count

## Details

This function runs a spatial accumulation curve starting from each site
individually, then extracts summary metrics from each curve. This allows
you to identify:

- Sites in species-rich areas (high initial slope)

- Core vs edge sites (fast vs slow accumulation)

- Spatial patterns in community structure

The metrics can be plotted as a heatmap using
`plot(result, type = "heatmap")`, which requires the `ggplot2` package.
For more sophisticated spatial visualization with study area boundaries,
see the `areaOfEffect` package.

## References

Soberon, J.M. & Llorente, J.B. (1993). The use of species accumulation
functions for the prediction of species richness. Conservation Biology,
7, 480-488.

## See also

[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) for
standard accumulation curves

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute per-site metrics
metrics <- spaccMetrics(species, coords,
                        metrics = c("slope_10", "half_richness", "auc"))

# Basic heatmap
plot(metrics, metric = "slope_10", type = "heatmap")

# Access metric values directly
metrics$metrics$slope_10
} # }
```
