#' Coverage-Based Spatial Rarefaction
#'
#' Compute spatial accumulation curves with sample coverage tracking.
#' Allows standardization by completeness (coverage) rather than sample size,
#' following Chao & Jost (2012).
#'
#' @param x A site-by-species matrix with abundance data.
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method. Default `"knn"`.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final coverage and richness for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#'
#' @return An object of class `spacc_coverage` containing:
#'   \item{richness}{Matrix of species richness (n_seeds x n_sites)}
#'   \item{individuals}{Matrix of individual counts}
#'   \item{coverage}{Matrix of coverage estimates}
#'   \item{coords, n_seeds, n_sites, method}{Parameters used}
#'
#' @details
#' Sample coverage (Chao & Jost 2012) estimates the proportion of the total
#' community abundance represented by observed species. It provides a measure
#' of sampling completeness that is independent of sample size.
#'
#' Coverage-based rarefaction allows fair comparison of diversity across
#' communities with different abundances by standardizing to equal completeness
#' rather than equal sample size.
#'
#' @references
#' Chao, A. & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. Ecology, 93, 2533-2547.
#'
#' @seealso [iNEXT::iNEXT()] for coverage-based rarefaction without spatial structure
#'
#' @examples
#' \dontrun{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' cov <- spaccCoverage(species, coords)
#' plot(cov)
#'
#' # Interpolate richness at 90% and 95% coverage
#' interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
#' }
#'
#' @export
spaccCoverage <- function(x,
                          coords,
                          n_seeds = 50L,
                          method = "knn",
                          distance = c("euclidean", "haversine"),
                          parallel = TRUE,
                          n_cores = NULL,
                          progress = TRUE,
                          seed = NULL,
                          map = FALSE) {

  distance <- match.arg(distance)

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    if (progress) cli_info("Computing site distances")
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)

  # Need abundance data for coverage calculation
  storage.mode(x) <- "integer"

  if (progress) cli_info(sprintf("Computing coverage-based accumulation (%d seeds)", n_seeds))

  result <- cpp_knn_coverage_parallel(x, dist_mat, n_seeds, n_cores, progress)

  if (progress) cli_success("Done")

  # Compute per-site map values if requested
  site_values <- NULL
  if (map) {
    if (progress) cli_info("Computing per-site coverage map values (all sites as seeds)")
    map_result <- cpp_knn_coverage_parallel(x, dist_mat, n_sites, n_cores, progress)

    site_values <- data.frame(
      site_id = seq_len(n_sites),
      x = coord_data$x,
      y = coord_data$y,
      final_richness = map_result$richness[, n_sites],
      final_coverage = map_result$coverage[, n_sites],
      final_individuals = map_result$individuals[, n_sites]
    )
    if (progress) cli_success("Map values computed")
  }

  structure(
    list(
      richness = result$richness,
      individuals = result$individuals,
      coverage = result$coverage,
      coords = coord_data,
      site_values = site_values,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_coverage"
  )
}


#' Interpolate Richness at Target Coverage Levels
#'
#' Estimate species richness at specified coverage levels by interpolation.
#'
#' @param x A `spacc_coverage` object.
#' @param target Numeric vector of target coverage levels (0 to 1).
#'   Default `c(0.90, 0.95, 0.99)`.
#'
#' @return A data.frame with columns for each target coverage level,
#'   showing interpolated richness for each seed.
#'
#' @export
interpolateCoverage <- function(x, target = c(0.90, 0.95, 0.99)) {
  stopifnot(inherits(x, "spacc_coverage"))
  stopifnot(all(target >= 0 & target <= 1))

  n_seeds <- x$n_seeds
  result <- matrix(NA, nrow = n_seeds, ncol = length(target))
  colnames(result) <- paste0("C", target * 100)

  for (s in 1:n_seeds) {
    richness <- x$richness[s, ]
    coverage <- x$coverage[s, ]

    result[s, ] <- interpolate_at_coverage(richness, coverage, target)
  }

  as.data.frame(result)
}


#' @export
print.spacc_coverage <- function(x, ...) {
  cat(sprintf("spacc coverage: %d sites, %d species, %d seeds\n",
              x$n_sites, x$n_species, x$n_seeds))

  mean_final_cov <- mean(x$coverage[, x$n_sites])
  mean_final_rich <- mean(x$richness[, x$n_sites])

  cat(sprintf("Mean final coverage: %.1f%%\n", mean_final_cov * 100))
  cat(sprintf("Mean final richness: %.1f species\n", mean_final_rich))
  invisible(x)
}


#' @export
summary.spacc_coverage <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2

  data.frame(
    sites = 1:object$n_sites,
    mean_richness = colMeans(object$richness),
    richness_lower = apply(object$richness, 2, stats::quantile, alpha),
    richness_upper = apply(object$richness, 2, stats::quantile, 1 - alpha),
    mean_individuals = colMeans(object$individuals),
    mean_coverage = colMeans(object$coverage),
    coverage_lower = apply(object$coverage, 2, stats::quantile, alpha),
    coverage_upper = apply(object$coverage, 2, stats::quantile, 1 - alpha)
  )
}


#' @export
plot.spacc_coverage <- function(x, type = c("curve", "map"),
                                 xaxis = c("sites", "coverage", "individuals"),
                                 ci = TRUE, ci_alpha = 0.2,
                                 metric = c("final_coverage", "final_richness", "final_individuals"),
                                 point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$site_values)) stop("No map data. Rerun spaccCoverage() with map = TRUE.")
    metric <- match.arg(metric)
    return(plot_spatial_map(x$site_values, metric,
                            title = sprintf("Coverage map: %s", metric),
                            subtitle = sprintf("%d sites, %s method", x$n_sites, x$method),
                            point_size = point_size, palette = palette))
  }

  check_suggests("ggplot2")

  xaxis <- match.arg(xaxis)
  summ <- summary(x)

  if (xaxis == "sites") {
    summ$x <- summ$sites
    xlab <- "Sites accumulated"
  } else if (xaxis == "coverage") {
    summ$x <- summ$mean_coverage
    xlab <- "Sample coverage"
  } else {
    summ$x <- summ$mean_individuals
    xlab <- "Individuals sampled"
  }

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = x, y = mean_richness))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = richness_lower, ymax = richness_upper),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::labs(
      x = xlab,
      y = "Species richness",
      title = "Coverage-Based Spatial Accumulation"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' @export
as_sf.spacc_coverage <- function(x, crs = NULL) {
  if (is.null(x$site_values)) stop("No map data. Rerun spaccCoverage() with map = TRUE.")
  as_sf_from_df(x$site_values, crs = crs)
}
