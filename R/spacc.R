#' Spatial Species Accumulation Curves
#'
#' Compute species accumulation curves using various spatial sampling methods
#' with C++ backend for performance.
#'
#' @param x A site-by-species matrix (rows = sites, cols = species) with
#'   presence/absence (0/1) or abundance data. Can also be a data.frame.
#' @param coords A data.frame with columns `x` and `y` containing site coordinates,
#'   or a `spacc_dist` object from [distances()].
#' @param n_seeds Integer. Number of random starting points for uncertainty
#'   quantification. Default 50.
#' @param method Character. Accumulation method:
#'   - `"knn"`: k-Nearest Neighbor (always visit closest unvisited)
#'   - `"kncn"`: k-Nearest Centroid Neighbor (visit closest to centroid)
#'   - `"random"`: Random order (null model)
#'   - `"radius"`: Expand by distance from seed
#'   - `"gaussian"`: Probabilistic selection weighted by distance
#'   - `"cone"`: Directional expansion within angular constraint
#'   - `"collector"`: Sites in data order (no randomization, single curve)
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param sigma Numeric. Bandwidth for Gaussian method. Default auto-calculated.
#' @param cone_width Numeric. Half-width in radians for cone method. Default pi/4.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores. Default `NULL` uses `detectCores() - 1`.
#' @param progress Logical. Show progress bar? Default `TRUE`.
#' @param seed Integer. Random seed for reproducibility. Default `NULL`.
#'
#' @return An object of class `spacc` containing:
#'   \item{curves}{Matrix of cumulative species counts (n_seeds x n_sites)}
#'   \item{coords}{Original coordinates}
#'   \item{n_seeds}{Number of seeds used}
#'   \item{method}{Method used}
#'   \item{n_species}{Total species in dataset}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' sac <- spacc(species, coords)
#' plot(sac)
#'
#' # Different methods
#' sac_knn <- spacc(species, coords, method = "knn")
#' sac_gauss <- spacc(species, coords, method = "gaussian", sigma = 10)
#' sac_cone <- spacc(species, coords, method = "cone", cone_width = pi/6)
#'
#' # Compare to null model
#' sac_rand <- spacc(species, coords, method = "random")
#' comp <- compare(sac_knn, sac_rand)
#' }
#'
#' @export
spacc <- function(x,
                  coords,
                  n_seeds = 50L,
                  method = c("knn", "kncn", "random", "radius", "gaussian", "cone", "collector"),
                  distance = c("euclidean", "haversine"),
                  sigma = NULL,
                  cone_width = pi / 4,
                  parallel = TRUE,
                  n_cores = NULL,
                  progress = TRUE,
                  seed = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)

  # Set RNG seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Resolve cores
  n_cores <- resolve_cores(n_cores, parallel)

  # Input validation
  x <- as.matrix(x)

  # Handle coords: either data.frame or spacc_dist
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot(
      "coords must have x and y columns" = all(c("x", "y") %in% names(coords))
    )
    coord_data <- coords
    dist_mat <- NULL
  }

  stopifnot(
    "x and coords must have same number of rows" = nrow(x) == nrow(coord_data),
    "n_seeds must be positive integer" = n_seeds > 0
  )

  n_sites <- nrow(x)
  n_species_total <- ncol(x)

  # Convert to presence/absence if abundance
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  # Collector method: no simulation needed
  if (method == "collector") {
    curve <- cpp_collector_single(species_pa)
    curves <- matrix(curve, nrow = 1)
    n_seeds <- 1L
  } else {
    # Compute distance matrix if not provided (needed for most methods)
    if (is.null(dist_mat) && method %in% c("knn", "radius", "gaussian")) {
      if (progress) cli_info(sprintf("Computing distances (%d x %d)", n_sites, n_sites))
      dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
    }

    # Auto-calculate sigma for Gaussian method
    if (method == "gaussian" && is.null(sigma)) {
      if (is.null(dist_mat)) {
        dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
      }
      sigma <- stats::median(dist_mat[dist_mat > 0])
    }

    # Run accumulation curves
    if (progress) cli_info(sprintf("Running %s accumulation (%d seeds, %d cores)", method, n_seeds, n_cores))

    curves <- switch(method,
      knn = cpp_knn_parallel(species_pa, dist_mat, n_seeds, n_cores, progress),
      kncn = cpp_kncn_parallel(species_pa, coord_data$x, coord_data$y, n_seeds, n_cores, progress),
      random = cpp_random_parallel(species_pa, n_seeds, n_cores, progress),
      radius = cpp_radius_parallel(species_pa, dist_mat, n_seeds, n_cores, progress),
      gaussian = cpp_gaussian_parallel(species_pa, dist_mat, n_seeds, sigma, n_cores, progress),
      cone = cpp_cone_parallel(species_pa, coord_data$x, coord_data$y, n_seeds, cone_width, n_cores, progress)
    )
  }

  if (progress) cli_success("Done")

  structure(
    list(
      curves = curves,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species_total,
      method = method,
      distance = distance,
      sigma = sigma,
      cone_width = if (method == "cone") cone_width else NULL,
      call = match.call()
    ),
    class = "spacc"
  )
}


#' Wavefront Expansion Accumulation
#'
#' Accumulate species within an expanding radius from seed points.
#' Models invasion spread from introduction points.
#'
#' @param x A site-by-species matrix.
#' @param coords A data.frame with x and y columns, or a `spacc_dist` object.
#' @param n_seeds Integer. Number of random starting points.
#' @param r0 Numeric. Initial radius. Default 0.
#' @param dr Numeric. Radius increment per step. Default auto-calculated.
#' @param n_steps Integer. Number of expansion steps. Default 50.
#' @param distance Character. Distance method.
#' @param progress Logical. Show progress?
#' @param seed Integer. Random seed.
#'
#' @return An object of class `spacc_wavefront` containing:
#'   \item{curves}{Matrix of species counts (n_seeds x n_steps)}
#'   \item{radius}{Vector of radius values}
#'   \item{sites_included}{Matrix of sites included at each step}
#'
#' @examples
#' \dontrun{
#' wf <- wavefront(species, coords, n_seeds = 20, n_steps = 100)
#' plot(wf)
#' }
#'
#' @export
wavefront <- function(x,
                      coords,
                      n_seeds = 50L,
                      r0 = 0,
                      dr = NULL,
                      n_steps = 50L,
                      distance = c("euclidean", "haversine"),
                      progress = TRUE,
                      seed = NULL) {

  distance <- match.arg(distance)

  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  n_sites <- nrow(x)
  n_species_total <- ncol(x)

  # Auto-calculate dr if not provided
  if (is.null(dr)) {
    max_dist <- max(dist_mat)
    dr <- max_dist / n_steps
  }

  # Convert to presence/absence
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  if (progress) cli_info(sprintf("Running wavefront expansion (%d seeds, %d steps)", n_seeds, n_steps))

  result <- cpp_wavefront_parallel(species_pa, dist_mat, n_seeds, r0, dr, n_steps, 1L, progress)

  if (progress) cli_success("Done")

  structure(
    list(
      curves = result$curves,
      radius = result$radius,
      sites_included = result$sites_included,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species_total,
      r0 = r0,
      dr = dr,
      n_steps = n_steps,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_wavefront"
  )
}


#' @export
print.spacc_wavefront <- function(x, ...) {
  cat(sprintf("spacc wavefront: %d sites, %d species, %d seeds\n",
              x$n_sites, x$n_species, x$n_seeds))
  cat(sprintf("Radius: %.2f to %.2f (%d steps)\n",
              x$r0, max(x$radius), x$n_steps))
  invisible(x)
}


#' @export
summary.spacc_wavefront <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2

  mean_curve <- colMeans(object$curves)
  mean_sites <- colMeans(object$sites_included)
  lower <- apply(object$curves, 2, stats::quantile, probs = alpha)
  upper <- apply(object$curves, 2, stats::quantile, probs = 1 - alpha)

  data.frame(
    radius = object$radius,
    mean_species = mean_curve,
    lower = lower,
    upper = upper,
    mean_sites = mean_sites
  )
}


#' @export
plot.spacc_wavefront <- function(x, ci = TRUE, ci_alpha = 0.3,
                                  xaxis = c("radius", "sites"), ...) {
  check_suggests("ggplot2")

  xaxis <- match.arg(xaxis)
  summ <- summary(x)

  if (xaxis == "radius") {
    df <- data.frame(
      x = summ$radius,
      mean = summ$mean_species,
      lower = summ$lower,
      upper = summ$upper
    )
    xlab <- "Radius"
  } else {
    df <- data.frame(
      x = summ$mean_sites,
      mean = summ$mean_species,
      lower = summ$lower,
      upper = summ$upper
    )
    xlab <- "Sites included"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = mean))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::labs(
      x = xlab,
      y = "Cumulative species",
      title = "Wavefront Species Accumulation",
      subtitle = sprintf("%d seeds, r0=%.1f, dr=%.2f", x$n_seeds, x$r0, x$dr)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Distance-Decay Analysis
#'
#' Analyze how species richness changes with distance from focal points.
#'
#' @param x A site-by-species matrix.
#' @param coords A data.frame with x and y columns.
#' @param n_seeds Integer. Number of focal points.
#' @param breaks Numeric vector. Distance thresholds. Default auto-calculated.
#' @param distance Character. Distance method.
#' @param progress Logical. Show progress?
#' @param seed Integer. Random seed.
#'
#' @return An object of class `spacc_decay` with distance-species relationships.
#'
#' @export
distance_decay <- function(x,
                           coords,
                           n_seeds = 50L,
                           breaks = NULL,
                           distance = c("euclidean", "haversine"),
                           progress = TRUE,
                           seed = NULL) {

  distance <- match.arg(distance)

  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)
  stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))

  n_sites <- nrow(x)

  # Compute distance matrix
  dist_mat <- cpp_distance_matrix(coords$x, coords$y, distance)

  # Auto-calculate breaks if not provided
  if (is.null(breaks)) {
    max_dist <- max(dist_mat)
    breaks <- seq(0, max_dist, length.out = 51)[-1]  # exclude 0
  }

  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  if (progress) cli_info(sprintf("Computing distance-decay (%d seeds, %d distances)", n_seeds, length(breaks)))

  curves <- cpp_distance_decay_parallel(species_pa, dist_mat, n_seeds, breaks, 1L, progress)

  if (progress) cli_success("Done")

  structure(
    list(
      curves = curves,
      breaks = breaks,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = ncol(x),
      distance = distance,
      call = match.call()
    ),
    class = "spacc_decay"
  )
}


#' @export
print.spacc_decay <- function(x, ...) {
  cat(sprintf("spacc distance-decay: %d seeds, %d distance breaks\n",
              x$n_seeds, length(x$breaks)))
  cat(sprintf("Distance range: %.2f to %.2f\n", min(x$breaks), max(x$breaks)))
  invisible(x)
}


#' @export
plot.spacc_decay <- function(x, ci = TRUE, ci_alpha = 0.3, ...) {
  check_suggests("ggplot2")

  mean_curve <- colMeans(x$curves)
  lower <- apply(x$curves, 2, stats::quantile, 0.025)
  upper <- apply(x$curves, 2, stats::quantile, 0.975)

  df <- data.frame(
    distance = x$breaks,
    mean = mean_curve,
    lower = lower,
    upper = upper
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = distance, y = mean))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::labs(
      x = "Distance from focal point",
      y = "Cumulative species",
      title = "Distance-Decay Relationship",
      subtitle = sprintf("%d focal points", x$n_seeds)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Resolve number of cores
#' @noRd
resolve_cores <- function(n_cores, parallel) {
  if (!parallel) return(1L)
  if (is.null(n_cores)) {
    max(1L, parallel::detectCores() - 1L)
  } else {
    as.integer(n_cores)
  }
}
