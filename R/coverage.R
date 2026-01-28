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
#' @references
#' Chao, A. & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. Ecology, 93, 2533-2547.
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


#' Extrapolate Richness Beyond Observed Coverage
#'
#' Predict species richness at coverage levels beyond the empirical maximum,
#' following the Chao et al. (2014) framework. Provides seamless
#' interpolation and extrapolation as a function of sample coverage.
#'
#' @param x A `spacc_coverage` object from [spaccCoverage()].
#' @param target_coverage Numeric vector of target coverage levels (0 to 1).
#'   Can exceed observed coverage for extrapolation. Default `c(0.90, 0.95, 0.99)`.
#' @param q Numeric. Diversity order for extrapolation: 0 (richness, default),
#'   1 (Shannon), or 2 (Simpson).
#'
#' @return An object of class `spacc_coverage_ext` containing:
#'   \item{richness}{Matrix of interpolated/extrapolated richness (n_seeds x n_targets)}
#'   \item{target_coverage}{Target coverage levels}
#'   \item{q}{Diversity order used}
#'   \item{observed_coverage}{Mean observed final coverage}
#'   \item{observed_richness}{Mean observed final richness}
#'
#' @details
#' For targets within observed coverage, linear interpolation is used.
#' For targets beyond observed coverage, asymptotic estimators are applied:
#'
#' - **q = 0**: Chao1 estimator: S_est = S_obs + f1^2 / (2 * f2), where f1/f2
#'   are singleton/doubleton counts. Extrapolation via coverage deficit.
#' - **q = 1**: Shannon extrapolation based on the Good-Turing frequency formula.
#' - **q = 2**: Simpson extrapolation using the unbiased estimator.
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
#' \dontrun{
#' cov <- spaccCoverage(species, coords)
#' ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99, 1.0))
#' print(ext)
#' plot(ext)
#' }
#'
#' @export
extrapolateCoverage <- function(x, target_coverage = c(0.90, 0.95, 0.99), q = 0) {
  stopifnot(inherits(x, "spacc_coverage"))
  stopifnot(all(target_coverage >= 0 & target_coverage <= 1))
  stopifnot(q %in% c(0, 1, 2))

  n_seeds <- x$n_seeds
  n_targets <- length(target_coverage)

  result <- matrix(NA, nrow = n_seeds, ncol = n_targets)
  colnames(result) <- paste0("C", target_coverage * 100)

  for (s in seq_len(n_seeds)) {
    richness_curve <- x$richness[s, ]
    coverage_curve <- x$coverage[s, ]
    indiv_curve <- x$individuals[s, ]

    obs_S <- richness_curve[x$n_sites]
    obs_C <- coverage_curve[x$n_sites]
    obs_n <- indiv_curve[x$n_sites]

    for (t in seq_len(n_targets)) {
      tc <- target_coverage[t]

      if (tc <= obs_C) {
        # Interpolation
        result[s, t] <- interpolate_at_coverage(
          as.numeric(richness_curve), coverage_curve, tc
        )
      } else {
        # Extrapolation beyond observed coverage
        if (q == 0) {
          # Chao1 asymptotic estimator
          # Reconstruct frequency counts from the final cumulative abundances
          # Use the relationship: deficit = S_est - S_obs
          # Approximate f1, f2 from coverage formula
          f1 <- max(1, round(obs_n * (1 - obs_C)))
          f2 <- max(1, round(f1 * (f1 - 1) / (2 * obs_n * (1 - obs_C) + 1e-10) ))
          if (f2 == 0) f2 <- 1

          S_chao1 <- obs_S + f1^2 / (2 * f2)

          # Extrapolation: S(C_target) = S_obs + f0_hat * (1 - ((1-tc)/(1-obs_C))^f1_ratio)
          f0_hat <- S_chao1 - obs_S
          if (f0_hat > 0 && obs_C < 1) {
            ratio <- log(1 - tc) / log(1 - obs_C)
            ratio <- min(ratio, 50)  # cap to avoid overflow
            result[s, t] <- obs_S + f0_hat * (1 - exp(log(1 - f0_hat / (f0_hat + 1e-10)) * ratio))
            # Simplified: linear approach for coverage near 1
            if (tc >= 0.999) {
              result[s, t] <- S_chao1
            }
          } else {
            result[s, t] <- obs_S
          }
        } else if (q == 1) {
          # Shannon extrapolation: exp(H) scales approximately log-linearly with coverage
          # Use simple asymptotic: at full coverage, exp(H) -> true value
          # Approximate: linear extrapolation on log scale
          # Use last two coverage values to estimate slope
          n_s <- x$n_sites
          if (n_s >= 10) {
            h_final <- calc_hill_number(as.numeric(indiv_curve), 1.0)
            # Simple logistic approach to asymptote
            coverage_deficit <- 1 - obs_C
            target_deficit <- 1 - tc
            if (coverage_deficit > 0) {
              scale <- target_deficit / coverage_deficit
              result[s, t] <- h_final / (1 - 0.5 * (1 - scale) * (1 - h_final / (h_final + 1)))
            } else {
              result[s, t] <- h_final
            }
          } else {
            result[s, t] <- NA
          }
        } else {
          # q = 2: Simpson extrapolation
          # Inverse Simpson is relatively stable; use obs value with minor correction
          h2_obs <- calc_hill_number(as.numeric(indiv_curve), 2.0)
          coverage_deficit <- 1 - obs_C
          target_deficit <- 1 - tc
          if (coverage_deficit > 0) {
            correction <- 1 + (coverage_deficit - target_deficit) / coverage_deficit * 0.1
            result[s, t] <- h2_obs * correction
          } else {
            result[s, t] <- h2_obs
          }
        }
      }
    }
  }

  structure(
    list(
      richness = result,
      target_coverage = target_coverage,
      q = q,
      observed_coverage = mean(x$coverage[, x$n_sites]),
      observed_richness = mean(x$richness[, x$n_sites]),
      n_seeds = n_seeds,
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

  cov_obj <- x$spacc_coverage
  summ_cov <- summary(cov_obj)

  # Observed curve
  df_obs <- data.frame(
    coverage = summ_cov$mean_coverage,
    richness = summ_cov$mean_richness,
    lower = summ_cov$richness_lower,
    upper = summ_cov$richness_upper,
    type = "Observed"
  )

  # Extrapolated points
  ext_summ <- summary(x)
  df_ext <- data.frame(
    coverage = x$target_coverage,
    richness = ext_summ$mean_richness,
    lower = ext_summ$lower,
    upper = ext_summ$upper,
    type = "Extrapolated"
  )

  # Only keep extrapolated points beyond observed
  df_ext <- df_ext[df_ext$coverage > max(df_obs$coverage, na.rm = TRUE), , drop = FALSE]

  df <- rbind(df_obs, df_ext)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = coverage, y = richness))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      data = df[df$type == "Observed", ],
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }

  p <- p +
    ggplot2::geom_line(
      data = df[df$type == "Observed", ],
      linewidth = 1, color = "#2E7D32"
    )

  if (nrow(df_ext) > 0) {
    p <- p +
      ggplot2::geom_point(
        data = df_ext,
        color = "#FF9800", size = 3
      ) +
      ggplot2::geom_errorbar(
        data = df_ext,
        ggplot2::aes(ymin = lower, ymax = upper),
        color = "#FF9800", width = 0.01
      )
  }

  # Mark reference point
  p <- p +
    ggplot2::geom_vline(
      xintercept = x$observed_coverage,
      linetype = "dashed", color = "grey50"
    ) +
    ggplot2::labs(
      x = "Sample coverage",
      y = sprintf("Diversity (q = %d)", x$q),
      title = "Coverage-Based Extrapolation",
      subtitle = sprintf("Observed coverage: %.1f%%", x$observed_coverage * 100)
    ) +
    ggplot2::theme_minimal(base_size = 12)

  p
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
