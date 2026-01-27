#' Spatial Beta Diversity Accumulation
#'
#' Analyze how beta diversity changes as sites are accumulated spatially.
#' Partitions beta diversity into turnover (species replacement) and
#' nestedness (species loss) components following Baselga (2010).
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method. Default `"knn"`.
#' @param index Character. Dissimilarity index: `"sorensen"` (default) or `"jaccard"`.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final beta values for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#'
#' @return An object of class `spacc_beta` containing:
#'   \item{beta_total}{Matrix of total beta diversity (n_seeds x n_sites-1)}
#'   \item{beta_turnover}{Matrix of turnover component}
#'   \item{beta_nestedness}{Matrix of nestedness component}
#'   \item{distance}{Matrix of cumulative distances}
#'   \item{n_seeds, n_sites, method, index}{Parameters used}
#'
#' @details
#' At each step of spatial accumulation, beta diversity is calculated between
#' the accumulated species pool and the newly added site. This reveals how
#' species composition changes as you expand spatially.
#'
#' **Interpretation:**
#' - High turnover: New sites contribute different species (replacement)
#' - High nestedness: New sites contribute subsets of existing species (loss)
#'
#' The sum of turnover and nestedness equals total beta diversity.
#'
#' @references
#' Baselga, A. (2010). Partitioning the turnover and nestedness components
#' of beta diversity. Global Ecology and Biogeography, 19, 134-143.
#'
#' @seealso [betapart::beta.pair()] for pairwise beta diversity
#'
#' @examples
#' \dontrun{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#'
#' beta <- spaccBeta(species, coords, n_seeds = 30)
#' plot(beta)
#'
#' # Compare Sorensen vs Jaccard
#' beta_jac <- spaccBeta(species, coords, index = "jaccard")
#' }
#'
#' @export
spaccBeta <- function(x,
                      coords,
                      n_seeds = 50L,
                      method = "knn",
                      index = c("sorensen", "jaccard"),
                      distance = c("euclidean", "haversine"),
                      parallel = TRUE,
                      n_cores = NULL,
                      progress = TRUE,
                      seed = NULL,
                      map = FALSE) {

  index <- match.arg(index)
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

  # Convert to presence/absence
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  use_jaccard <- index == "jaccard"

  if (progress) cli_info(sprintf("Computing beta diversity (%s, %d seeds)",
                                  index, n_seeds))

  result <- cpp_beta_knn_parallel(species_pa, dist_mat, n_seeds,
                                   use_jaccard, n_cores, progress)

  if (progress) cli_success("Done")

  # Compute per-site map values if requested
  site_values <- NULL
  if (map) {
    if (progress) cli_info("Computing per-site beta map values (all sites as seeds)")
    map_result <- cpp_beta_knn_parallel(species_pa, dist_mat, n_sites,
                                         use_jaccard, n_cores, progress)

    n_steps <- ncol(map_result$beta_total)
    site_values <- data.frame(
      site_id = seq_len(n_sites),
      x = coord_data$x,
      y = coord_data$y,
      beta_total = map_result$beta_total[, n_steps],
      beta_turnover = map_result$beta_turnover[, n_steps],
      beta_nestedness = map_result$beta_nestedness[, n_steps]
    )
    if (progress) cli_success("Map values computed")
  }

  structure(
    list(
      beta_total = result$beta_total,
      beta_turnover = result$beta_turnover,
      beta_nestedness = result$beta_nestedness,
      distance = result$distance,
      coords = coord_data,
      site_values = site_values,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      index = index,
      call = match.call()
    ),
    class = "spacc_beta"
  )
}


#' @export
print.spacc_beta <- function(x, ...) {
  cat(sprintf("spacc beta diversity: %d sites, %d seeds\n",
              x$n_sites, x$n_seeds))
  cat(sprintf("Index: %s, Method: %s\n", x$index, x$method))

  # Final values
  final_total <- mean(x$beta_total[, ncol(x$beta_total)])
  final_turn <- mean(x$beta_turnover[, ncol(x$beta_turnover)])
  final_nest <- mean(x$beta_nestedness[, ncol(x$beta_nestedness)])

 cat(sprintf("Mean final beta: %.3f (turnover: %.3f, nestedness: %.3f)\n",
              final_total, final_turn, final_nest))
  invisible(x)
}


#' @export
summary.spacc_beta <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2
  n_steps <- ncol(object$beta_total)

  data.frame(
    step = 1:n_steps,
    mean_distance = colMeans(object$distance),
    beta_total = colMeans(object$beta_total),
    beta_total_lower = apply(object$beta_total, 2, stats::quantile, alpha),
    beta_total_upper = apply(object$beta_total, 2, stats::quantile, 1 - alpha),
    beta_turnover = colMeans(object$beta_turnover),
    beta_nestedness = colMeans(object$beta_nestedness)
  )
}


#' @export
plot.spacc_beta <- function(x, type = c("curve", "map"), partition = TRUE,
                            xaxis = c("sites", "distance"),
                            ci = TRUE, ci_alpha = 0.2,
                            component = c("beta_total", "beta_turnover", "beta_nestedness"),
                            point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$site_values)) stop("No map data. Rerun spaccBeta() with map = TRUE.")
    component <- match.arg(component)
    return(plot_spatial_map(x$site_values, component,
                            title = sprintf("Beta diversity map: %s", component),
                            subtitle = sprintf("%d sites, %s (%s)", x$n_sites, x$index, x$method),
                            point_size = point_size, palette = palette))
  }

  check_suggests("ggplot2")

  xaxis <- match.arg(xaxis)
  summ <- summary(x)

  if (xaxis == "sites") {
    summ$x <- summ$step + 1  # step 1 = site 2
    xlab <- "Sites accumulated"
  } else {
    summ$x <- summ$mean_distance
    xlab <- "Cumulative distance"
  }

  if (partition) {
    # Reshape for stacked/faceted plot
    df <- data.frame(
      x = rep(summ$x, 3),
      value = c(summ$beta_total, summ$beta_turnover, summ$beta_nestedness),
      component = rep(c("Total", "Turnover", "Nestedness"), each = nrow(summ))
    )
    df$component <- factor(df$component, levels = c("Total", "Turnover", "Nestedness"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, color = component)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_color_manual(values = c("Total" = "#2E7D32",
                                              "Turnover" = "#1565C0",
                                              "Nestedness" = "#C62828"))
  } else {
    p <- ggplot2::ggplot(summ, ggplot2::aes(x = x, y = beta_total))

    if (ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = beta_total_lower, ymax = beta_total_upper),
        alpha = ci_alpha, fill = "#4CAF50"
      )
    }

    p <- p + ggplot2::geom_line(linewidth = 1, color = "#2E7D32")
  }

  p +
    ggplot2::labs(
      x = xlab,
      y = paste0("Beta diversity (", x$index, ")
"),
      color = "Component",
      title = "Spatial Beta Diversity Accumulation"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}


#' @export
as_sf.spacc_beta <- function(x, crs = NULL) {
  if (is.null(x$site_values)) stop("No map data. Rerun spaccBeta() with map = TRUE.")
  as_sf_from_df(x$site_values, crs = crs)
}
