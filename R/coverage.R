#' Coverage-Based Spatial Rarefaction
#'
#' Compute spatial accumulation curves with sample coverage tracking.
#' Allows standardization by completeness (coverage) rather than sample size,
#' following Chao & Jost (2012) or the sample-based estimator of Chiu (2023).
#'
#' @param x A site-by-species matrix with abundance data.
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method: canonical fixed-focus `"knn"`
#'   or `"nn_walk"`. Default `"knn"`.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param coverage Character. Coverage estimator to use: `"chao"` (default)
#'   for the individual-based Chao & Jost (2012) estimator using
#'   singletons/doubletons, or `"chiu"` for the sample-based Chiu (2023)
#'   estimator using incidence frequency counts (Q1/Q2) and unique-species
#'   abundance (G1). The Chiu estimator is recommended for spatially aggregated
#'   data where sampling units are plots/sites rather than independent individuals.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final coverage and richness for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#'
#' @return An object of class `spacc_coverage` containing:
#'   \item{richness}{Matrix of species richness (n_seeds x n_sites)}
#'   \item{individuals}{Matrix of individual counts}
#'   \item{coverage}{Matrix of coverage estimates}
#'   \item{pooled}{Species abundances pooled over all sites}
#'   \item{coverage_type}{Coverage estimator used (`"chao"` or `"chiu"`)}
#'   \item{coords, n_seeds, n_sites, method}{Parameters used}
#'
#' @details
#' Sample coverage estimates the proportion of the total community abundance
#' represented by observed species. It provides a measure of sampling
#' completeness that is independent of sample size.
#'
#' The **Chao-Jost (2012)** estimator counts singletons (f1) and doubletons
#' (f2) in the cumulative abundance vector. It assumes individuals are sampled
#' independently, which may not hold for plot-based spatial data.
#'
#' The **Chiu (2023)** estimator uses incidence frequency counts instead:
#' Q1 (species in exactly 1 site), Q2 (species in exactly 2 sites), and G1
#' (total abundance of Q1 species). It gives near-unbiased coverage estimates
#' when organisms are spatially aggregated across sampling units.
#'
#' Coverage-based rarefaction allows fair comparison of diversity across
#' communities with different abundances by standardizing to equal completeness
#' rather than equal sample size.
#'
#' @references
#' Chao, A. & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. Ecology, 93, 2533-2547.
#'
#' Chiu, C.-H. (2023). A sample-based estimator for sample coverage.
#' Ecology, 104, e4099.
#'
#' @seealso [iNEXT::iNEXT()] for coverage-based rarefaction without spatial structure
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' cov <- spaccCoverage(species, coords)
#' plot(cov)
#'
#' # Sample-based coverage (recommended for spatial data)
#' cov_chiu <- spaccCoverage(species, coords, coverage = "chiu")
#'
#' # Interpolate richness at 90% and 95% coverage
#' interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
#' }
#'
#' @export
spaccCoverage <- function(x,
                          coords,
                          n_seeds = 50L,
                          method = c("knn", "nn_walk"),
                          distance = c("euclidean", "haversine"),
                          coverage = c("chao", "chiu"),
                          parallel = TRUE,
                          n_cores = NULL,
                          progress = TRUE,
                          seed = NULL,
                          map = FALSE,
                          focal_points = NULL,
                          focal_domain = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  coverage <- match.arg(coverage)

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  site <- .resolve_coords(coords)
  coord_data <- site$coords
  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)
  storage.mode(x) <- "integer"

  coverage_type <- if (coverage == "chiu") 1L else 0L
  cov_label <- if (coverage == "chiu") "Chiu 2023" else "Chao-Jost 2012"
  if (progress) cli_info(sprintf("Computing coverage-based accumulation (%d seeds, %s)", n_seeds, cov_label))

  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   site$dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)
  run <- function(orders) {
    cpp_order_coverage_parallel(x, orders - 1L, coverage_type, n_cores)
  }
  result <- run(ordering$orders)

  if (progress) cli_success("Done")

  site_values <- NULL
  if (map) {
    if (progress) cli_info("Computing per-site coverage map values (all sites as seeds)")
    site_values <- .site_map_values(method, coord_data, distance, site$dist_mat, n_cores,
      function(orders) {
        map_result <- run(orders)
        list(final_richness = .final_per_seed(map_result$richness),
             final_coverage = .final_per_seed(map_result$coverage),
             final_individuals = .final_per_seed(map_result$individuals))
      })
    if (progress) cli_success("Map values computed")
  }

  structure(
    list(
      richness = result$richness,
      individuals = result$individuals,
      coverage = result$coverage,
      pooled = colSums(x),
      coords = coord_data,
      site_values = site_values,
      n_seeds = ordering$n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      distance = distance,
      coverage_type = coverage,
      focal_points = ordering$focal_points,
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
#' @references
#' Chao, A. & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. Ecology, 93, 2533-2547.
#'
#' @export
interpolateCoverage <- function(x, target = c(0.90, 0.95, 0.99)) {
  stopifnot(inherits(x, "spacc_coverage"))
  stopifnot(all(target >= 0 & target <= 1))

  rows <- lapply(seq_len(x$n_seeds), function(s) {
    interpolate_at_coverage(x$richness[s, ], x$coverage[s, ], target)
  })
  result <- matrix(unlist(rows), nrow = x$n_seeds, byrow = TRUE)
  colnames(result) <- paste0("C", target * 100)
  as.data.frame(result)
}


#' Extrapolate Diversity Beyond Observed Coverage
#'
#' Hill numbers at target sample coverages along a coverage-tracking
#' accumulation, interpolated within the observed coverage and extrapolated
#' beyond it (Chao et al. 2014).
#'
#' @param x A `spacc_coverage` object from [spaccCoverage()] (richness,
#'   `q = 0` only) or a `spacc_hill_coverage` object from
#'   [spaccHillCoverage()] (any of its orders).
#' @param target_coverage Numeric vector of target coverage levels (0 to 1).
#'   Can exceed observed coverage for extrapolation. Default `c(0.90, 0.95, 0.99)`.
#' @param q Numeric. Diversity order: 0 (richness, default), or for a
#'   `spacc_hill_coverage` object one of its orders.
#'
#' @return An object of class `spacc_coverage_ext` containing:
#'   \item{richness}{Matrix of interpolated/extrapolated Hill numbers of
#'     order `q` (n_seeds x n_targets)}
#'   \item{target_coverage}{Target coverage levels}
#'   \item{q}{Diversity order used}
#'   \item{observed_coverage}{Mean observed final coverage}
#'   \item{observed_richness}{Mean observed final Hill number of order `q`}
#'
#' @details
#' A target at or below a seed's final coverage is read off that seed's
#' spatial accumulation curve by linear interpolation in coverage. Above it,
#' every seed has accumulated the same pooled sample, so the value is the
#' coverage-based extrapolation of the pooled abundances (Chao et al. 2014):
#' the sample size reaching the target coverage is found from the
#' Chao-Jost coverage deficit, and the Hill number at that size approaches
#' the Chao1 (q = 0), Chao et al. (2013) Shannon (q = 1) or unbiased Simpson
#' (q = 2) asymptote; `target_coverage = 1` gives the asymptote itself.
#' Extrapolation is defined for integer `q` and for the individual-based
#' coverage (`coverage = "chao"`).
#'
#' @references
#' Chao, A. & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. Ecology, 93, 2533-2547.
#'
#' Chao, A., Gotelli, N.J., Hsieh, T.C., et al. (2014). Rarefaction and
#' extrapolation with Hill numbers: a framework for sampling and estimation in
#' species diversity studies. Ecological Monographs, 84, 45-67.
#'
#' @seealso [spaccCoverage()], [interpolateCoverage()]
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' cov <- spaccCoverage(species, coords)
#' ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99))
#' print(ext)
#' }
#'
#' @export
extrapolateCoverage <- function(x, target_coverage = c(0.90, 0.95, 0.99), q = 0) {
  stopifnot("x must be a spacc_coverage or spacc_hill_coverage object" =
              inherits(x, c("spacc_coverage", "spacc_hill_coverage")))
  stopifnot(all(target_coverage >= 0 & target_coverage <= 1))
  stopifnot("q must be a single number" = is.numeric(q) && length(q) == 1)

  curves <- if (inherits(x, "spacc_hill_coverage")) {
    if (!q %in% x$q) stop(sprintf("q = %g is not among the orders of x (%s)",
                                  q, paste(x$q, collapse = ", ")), call. = FALSE)
    x$hills[[match(q, x$q)]]
  } else {
    if (q != 0) stop("a spacc_coverage object carries richness only (q = 0); ",
                     "use spaccHillCoverage() for other orders", call. = FALSE)
    x$richness
  }
  if (is.null(x$pooled)) stop("x was created by an older spacc version; recompute it",
                              call. = FALSE)

  final_C <- x$coverage[, x$n_sites]
  beyond <- outer(final_C, target_coverage, `<`)
  if (any(beyond)) {
    if (identical(x$coverage_type, "chiu")) {
      stop("Extrapolation beyond the observed coverage uses the individual-based ",
           "estimators; recompute with coverage = \"chao\"", call. = FALSE)
    }
    if (q != round(q)) stop("Extrapolation is defined for integer q", call. = FALSE)
  }
  pooled <- calc_hill_at_coverage(as.numeric(x$pooled), target_coverage, q)$hill[1, ]

  result <- matrix(NA_real_, nrow = x$n_seeds, ncol = length(target_coverage))
  colnames(result) <- paste0("C", target_coverage * 100)
  for (s in seq_len(x$n_seeds)) {
    within <- !beyond[s, ]
    if (any(within)) {
      result[s, within] <- interpolate_at_coverage(as.numeric(curves[s, ]), x$coverage[s, ],
                                                   target_coverage[within])
    }
    result[s, !within] <- pooled[!within]
  }

  structure(
    list(
      richness = result,
      target_coverage = target_coverage,
      q = q,
      observed_coverage = mean(final_C),
      observed_richness = mean(curves[, x$n_sites]),
      n_seeds = x$n_seeds,
      spacc_coverage = x
    ),
    class = "spacc_coverage_ext"
  )
}


#' @export
print.spacc_coverage_ext <- function(x, ...) {
  cat("Coverage-based extrapolation\n")
  cat(strrep("-", 32), "\n")
  cat(sprintf("Diversity order: q = %d\n", x$q))
  cat(sprintf("Observed coverage: %.1f%%\n", x$observed_coverage * 100))
  cat(sprintf("Observed richness: %.1f\n", x$observed_richness))
  cat("\nExtrapolated richness:\n")
  means <- colMeans(x$richness, na.rm = TRUE)
  sds <- apply(x$richness, 2, stats::sd, na.rm = TRUE)
  for (i in seq_along(x$target_coverage)) {
    cat(sprintf("  C=%.0f%%: %.1f (+/- %.1f)\n",
                x$target_coverage[i] * 100, means[i], sds[i]))
  }
  invisible(x)
}


#' @export
summary.spacc_coverage_ext <- function(object, ...) {
  data.frame(
    target_coverage = object$target_coverage,
    mean_richness = colMeans(object$richness, na.rm = TRUE),
    sd = apply(object$richness, 2, stats::sd, na.rm = TRUE),
    lower = apply(object$richness, 2, stats::quantile, 0.025, na.rm = TRUE),
    upper = apply(object$richness, 2, stats::quantile, 0.975, na.rm = TRUE)
  )
}


#' @export
plot.spacc_coverage_ext <- function(x, ci = TRUE, ci_alpha = 0.2, ...) {
  check_suggests("ggplot2")

  summ_cov <- summary(x$spacc_coverage)
  df_obs <- data.frame(
    coverage = summ_cov$mean_coverage,
    richness = summ_cov$mean_richness,
    lower = summ_cov$richness_lower,
    upper = summ_cov$richness_upper
  )

  ext_summ <- summary(x)
  df_ext <- data.frame(
    coverage = x$target_coverage,
    richness = ext_summ$mean_richness,
    lower = ext_summ$lower,
    upper = ext_summ$upper
  )
  df_ext <- df_ext[df_ext$coverage > max(df_obs$coverage, na.rm = TRUE), , drop = FALSE]

  p <- .plot_curve(df_obs, x = "coverage", y = "richness",
                   xlab = "Sample coverage", ylab = sprintf("Diversity (q = %d)", x$q),
                   title = "Coverage-Based Extrapolation",
                   subtitle = sprintf("Observed coverage: %.1f%%", x$observed_coverage * 100),
                   ci = ci, ci_alpha = ci_alpha)

  if (nrow(df_ext) > 0) {
    p <- p +
      ggplot2::geom_point(data = df_ext, color = "#FF9800", size = 3) +
      ggplot2::geom_errorbar(data = df_ext,
                             ggplot2::aes(ymin = lower, ymax = upper),
                             color = "#FF9800", width = 0.01)
  }

  p + ggplot2::geom_vline(xintercept = x$observed_coverage,
                          linetype = "dashed", color = "grey50")
}


#' @export
print.spacc_coverage <- function(x, ...) {
  cov_label <- if (identical(x$coverage_type, "chiu")) "Chiu 2023" else "Chao-Jost 2012"
  .print_curves(x, sprintf("spacc coverage (%s)", cov_label),
                sprintf("Mean final coverage: %.1f%%", mean(x$coverage[, x$n_sites]) * 100),
                sprintf("Mean final richness: %.1f species", mean(x$richness[, x$n_sites])))
}


#' @export
summary.spacc_coverage <- function(object, ci_level = 0.95, ...) {
  richness <- .curve_stats(object$richness, ci_level)
  coverage <- .curve_stats(object$coverage, ci_level)
  data.frame(
    sites = seq_len(object$n_sites),
    mean_richness = richness$mean,
    richness_lower = richness$lower,
    richness_upper = richness$upper,
    mean_individuals = colMeans(object$individuals),
    mean_coverage = coverage$mean,
    coverage_lower = coverage$lower,
    coverage_upper = coverage$upper
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
    metric <- match.arg(metric)
    return(.plot_site_map(x, metric,
                          title = sprintf("Coverage map: %s", metric),
                          subtitle = sprintf("%d sites, %s method", x$n_sites, x$method),
                          point_size = point_size, palette = palette))
  }

  xaxis <- match.arg(xaxis)
  axis <- switch(xaxis,
    sites = c("sites", "Sites accumulated"),
    coverage = c("mean_coverage", "Sample coverage"),
    individuals = c("mean_individuals", "Individuals sampled")
  )
  .plot_curve(summary(x), x = axis[1], y = "mean_richness",
              lower = "richness_lower", upper = "richness_upper",
              xlab = axis[2], ylab = "Species richness",
              title = "Coverage-Based Spatial Accumulation",
              ci = ci, ci_alpha = ci_alpha)
}


#' @export
as_sf.spacc_coverage <- function(x, crs = NULL) .as_sf_site_values(x, crs)
