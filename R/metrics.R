#' Per-Site Accumulation Metrics
#'
#' Compute spatial accumulation metrics for each site as a starting point.
#' Useful for identifying sites with high or low accumulation rates,
#' visualizing spatial patterns in diversity, and understanding edge effects.
#'
#' @param x A site-by-species matrix (rows = sites, cols = species).
#' @param coords A data.frame with columns `x` and `y` containing site coordinates,
#'   or a `spacc_dist` object from [distances()].
#' @param metrics Character vector. Metrics to compute. Options include:
#'   `"slope_10"` (initial slope, first 10% of sites),
#'   `"slope_25"` (initial slope, first 25% of sites),
#'   `"half_richness"` (sites to reach 50% of total species),
#'   `"richness_50pct"` (alias for half_richness),
#'   `"richness_75pct"` (sites to reach 75% of species),
#'   `"richness_90pct"` (sites to reach 90% of species),
#'   `"auc"` (area under accumulation curve),
#'   `"final_richness"` (total species starting from this site).
#' @param method Character. Accumulation method: `"knn"`, `"kncn"`,
#'   `"nn_walk"`, or `"random"`.
#'   Default `"knn"`.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores for parallel processing.
#' @param progress Logical. Show progress bar? Default `TRUE`.
#'
#' @return An object of class `spacc_metrics` containing:
#'   \item{metrics}{Data frame with one row per site and columns for each metric}
#'   \item{coords}{Original coordinates}
#'   \item{metric_names}{Names of computed metrics}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Total species count}
#'
#' @details
#' This function runs a spatial accumulation curve starting from each site
#' individually, then extracts summary metrics from each curve. This allows
#' you to identify:
#' - Sites in species-rich areas (high initial slope)
#' - Core vs edge sites (fast vs slow accumulation)
#' - Spatial patterns in community structure
#'
#' The metrics can be plotted as a heatmap using `plot(result, type = "heatmap")`,
#' which requires the `ggplot2` package. For more sophisticated spatial
#' visualization with study area boundaries, see the `areaOfEffect` package.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' metrics <- spaccMetrics(species, coords,
#'                         metrics = c("slope_10", "auc"))
#' metrics$metrics$slope_10
#' }
#'
#' @references
#' Soberon, J.M. & Llorente, J.B. (1993). The use of species accumulation
#' functions for the prediction of species richness. Conservation Biology,
#' 7, 480-488.
#'
#' @seealso [spacc()] for standard accumulation curves
#' @export
spaccMetrics <- function(x,
                         coords,
                         metrics = c("slope_10", "half_richness", "auc"),
                         method = c("knn", "kncn", "nn_walk", "random"),
                         distance = c("euclidean", "haversine"),
                         parallel = TRUE,
                         n_cores = NULL,
                         progress = TRUE) {

  method <- match.arg(method)
  distance <- match.arg(distance)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  site <- .resolve_coords(coords)
  coord_data <- site$coords
  stopifnot("x and coords must have same number of rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)

  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  if (progress) cli_info(sprintf("Running %s accumulation from each site (%d sites, %d cores)",
                                  method, n_sites, n_cores))

  ordering <- .accumulation_orders(method, coord_data, n_sites, distance,
                                   site$dist_mat, all_sites = TRUE, n_cores = n_cores)
  curves <- cpp_order_parallel(species_pa, ordering$orders - 1L, n_cores)

  if (progress) cli_info("Extracting metrics")

  metric_df <- extract_metrics(curves, metrics, n_species)
  metric_df$site_id <- seq_len(n_sites)
  metric_df$x <- coord_data$x
  metric_df$y <- coord_data$y

  if (progress) cli_success("Done")

  structure(
    list(
      metrics = metric_df,
      coords = coord_data,
      curves = curves,
      metric_names = metrics,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_metrics"
  )
}


#' Extract metrics from accumulation curves
#' @noRd
extract_metrics <- function(curves, metrics, n_species) {
  n_sites <- nrow(curves)
  result <- data.frame(row.names = seq_len(n_sites))

  for (m in metrics) {
    result[[m]] <- switch(m,
      "slope_10" = calc_slope(curves, frac = 0.10),
      "slope_25" = calc_slope(curves, frac = 0.25),
      "half_richness" = calc_sites_to_richness(curves, frac = 0.50, n_species),
      "richness_50pct" = calc_sites_to_richness(curves, frac = 0.50, n_species),
      "richness_75pct" = calc_sites_to_richness(curves, frac = 0.75, n_species),
      "richness_90pct" = calc_sites_to_richness(curves, frac = 0.90, n_species),
      "auc" = calc_auc(curves),
      "final_richness" = curves[, ncol(curves)],
      stop("Unknown metric: ", m)
    )
  }

  result
}


#' Calculate initial slope
#' @noRd
calc_slope <- function(curves, frac) {
  n_sites <- ncol(curves)
  n_points <- max(2, round(n_sites * frac))

  apply(curves, 1, function(curve) {
    x <- seq_len(n_points)
    y <- curve[seq_len(n_points)]
    if (length(unique(y)) == 1) return(0)
    coef(lm(y ~ x))[2]
  })
}


#' Calculate sites to reach richness fraction
#' @noRd
calc_sites_to_richness <- function(curves, frac, n_species) {
  target <- n_species * frac

  apply(curves, 1, function(curve) {
    idx <- which(curve >= target)
    if (length(idx) == 0) return(length(curve))
    min(idx)
  })
}


#' Calculate area under curve
#' @noRd
calc_auc <- function(curves) {
  apply(curves, 1, function(curve) {
    sum(curve) / length(curve)
  })
}


#' @export
print.spacc_metrics <- function(x, ...) {
  cat(sprintf("spacc_metrics: %d sites, %d species\n", x$n_sites, x$n_species))
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Metrics: %s\n", paste(x$metric_names, collapse = ", ")))
  invisible(x)
}


#' @export
summary.spacc_metrics <- function(object, ...) {
  cat("Metric summary:\n")
  for (m in object$metric_names) {
    vals <- object$metrics[[m]]
    cat(sprintf("  %s: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
                m, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
                min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
  }
  invisible(object)
}


#' Internal spatial heatmap plotting helper
#' @noRd
plot_spatial_map <- function(df, value_col, title, subtitle = NULL,
                             point_size = 3, palette = "viridis") {
  check_suggests("ggplot2")

  ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], color = .data[[value_col]])) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_viridis_c(option = palette) +
    ggplot2::labs(
      x = "X coordinate",
      y = "Y coordinate",
      color = value_col,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::coord_equal() +
    spacc_theme() +
    ggplot2::theme(legend.position = "right")
}


#' Internal as_sf conversion helper
#' @noRd
as_sf_from_df <- function(df, crs = NULL) {
  check_suggests("sf")
  if (!is.null(crs)) {
    sf::st_as_sf(df, coords = c("x", "y"), crs = crs)
  } else {
    sf::st_as_sf(df, coords = c("x", "y"))
  }
}


#' Plot spacc_metrics
#'
#' Create visualizations of per-site accumulation metrics.
#'
#' @param x A `spacc_metrics` object from [spaccMetrics()].
#' @param metric Character. Which metric to plot. Default is first metric.
#' @param type Character. Plot type:
#'   \describe{
#'     \item{`"heatmap"`}{Spatial heatmap colored by metric value}
#'     \item{`"points"`}{Simple point plot (same as heatmap but clearer name)}
#'     \item{`"histogram"`}{Distribution of metric values}
#'   }
#' @param point_size Numeric. Size of points in heatmap. Default 3.
#' @param palette Character. Color palette for heatmap. One of `"magma"` (A),
#'   `"inferno"` (B), `"plasma"` (C), `"viridis"` (D, default), `"cividis"` (E),
#'   `"rocket"` (F), `"mako"` (G), `"turbo"` (H). Single letters also accepted.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' metrics <- spaccMetrics(species, coords, metrics = c("slope_10", "auc"))
#' plot(metrics, metric = "slope_10", type = "heatmap")
#' }
#'
#' @export
plot.spacc_metrics <- function(x, metric = NULL, type = c("heatmap", "points", "histogram"),
                                point_size = 3, palette = "viridis", ...) {

  check_suggests("ggplot2")

  type <- match.arg(type)

  if (is.null(metric)) {
    metric <- x$metric_names[1]
  }

  if (!metric %in% names(x$metrics)) {
    stop("Metric '", metric, "' not found. Available: ",
         paste(x$metric_names, collapse = ", "))
  }

  df <- x$metrics
  subtitle <- sprintf("%d sites, %s method", x$n_sites, x$method)

  if (type == "histogram") {
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data[[metric]])) +
        ggplot2::geom_histogram(bins = 30, fill = "#4CAF50", color = "white", alpha = 0.8) +
        ggplot2::labs(x = metric, y = "Count",
                      title = sprintf("Distribution: %s", metric), subtitle = subtitle) +
        spacc_theme()
    )
  }

  plot_spatial_map(df, metric,
                   title = sprintf("Spatial pattern: %s", metric),
                   subtitle = subtitle,
                   point_size = point_size, palette = palette)
}


#' Convert spacc_metrics to sf
#'
#' Convert metrics to an sf object for spatial analysis and integration
#' with the areaOfEffect package.
#'
#' @param x A `spacc_metrics` object.
#' @param crs Coordinate reference system. Default `NULL` (no CRS).
#'   Use EPSG codes like `4326` for WGS84 or `32631` for UTM zone 31N.
#'
#' @return An sf object with POINT geometries and metric columns.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' metrics <- spaccMetrics(species, coords)
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   metrics_sf <- as_sf(metrics)
#' }
#' }
#'
#' @export
as_sf <- function(x, crs = NULL) {
  UseMethod("as_sf")
}


#' @export
as_sf.spacc_metrics <- function(x, crs = NULL) {
  as_sf_from_df(x$metrics, crs = crs)
}
