# Plot spacc_metrics

Create visualizations of per-site accumulation metrics.

## Usage

``` r
# S3 method for class 'spacc_metrics'
plot(
  x,
  metric = NULL,
  type = c("heatmap", "points", "histogram"),
  point_size = 3,
  palette = "viridis",
  ...
)
```

## Arguments

- x:

  A `spacc_metrics` object from
  [`spaccMetrics()`](https://gillescolling.com/spacc/reference/spaccMetrics.md).

- metric:

  Character. Which metric to plot. Default is first metric.

- type:

  Character. Plot type:

  `"heatmap"`

  :   Spatial heatmap colored by metric value

  `"points"`

  :   Simple point plot (same as heatmap but clearer name)

  `"histogram"`

  :   Distribution of metric values

- point_size:

  Numeric. Size of points in heatmap. Default 3.

- palette:

  Character. Color palette for heatmap. Default `"viridis"`.

- ...:

  Additional arguments (unused).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
metrics <- spaccMetrics(species, coords, metrics = c("slope_10", "auc"))

# Heatmap of initial accumulation rate
plot(metrics, metric = "slope_10", type = "heatmap")

# Distribution of AUC values
plot(metrics, metric = "auc", type = "histogram")
} # }
```
