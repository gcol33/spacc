#' Spatial Phylogenetic Diversity Accumulation
#'
#' Compute spatial accumulation of phylogenetic diversity metrics (MPD, MNTD, PD).
#'
#' @param x A site-by-species matrix.
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param tree A phylogenetic tree of class `phylo` (from ape package), or
#'   a pairwise phylogenetic distance matrix.
#' @param metric Character vector. Metrics to compute:
#'   - `"mpd"`: Mean Pairwise Distance
#'   - `"mntd"`: Mean Nearest Taxon Distance
#'   - `"pd"`: Faith's Phylogenetic Diversity (requires tree, not distance matrix)
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method. Default `"knn"`.
#' @param distance Character. Site distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final values for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#'
#' @return An object of class `spacc_phylo` containing:
#'   \item{curves}{Named list of matrices, one per metric (n_seeds x n_sites)}
#'   \item{metric}{Metrics computed}
#'   \item{coords, n_seeds, n_sites, method}{Parameters used}
#'
#' @details
#' Phylogenetic diversity metrics incorporate evolutionary relationships:
#'
#' - **MPD (Mean Pairwise Distance)**: Average phylogenetic distance between
#'   all pairs of species. Sensitive to tree-wide patterns.
#'
#' - **MNTD (Mean Nearest Taxon Distance)**: Average distance to closest
#'   relative. Sensitive to terminal clustering.
#'
#' - **PD (Faith's Phylogenetic Diversity)**: Total branch length connecting
#'   species. Requires full tree object.
#'
#' @references
#' Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity.
#' Biological Conservation, 61, 1-10.
#'
#' Webb, C.O. (2000). Exploring the phylogenetic structure of ecological
#' communities: an example for rain forest trees. American Naturalist, 156, 145-155.
#'
#' @seealso [picante::mpd()], [picante::mntd()], [picante::pd()]
#'
#' @examples
#' \dontrun{
#' library(ape)
#'
#' # Create random tree
#' tree <- rtree(30)
#'
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' colnames(species) <- tree$tip.label
#'
#' phylo <- spaccPhylo(species, coords, tree, metric = c("mpd", "mntd"))
#' plot(phylo)
#' }
#'
#' @export
spaccPhylo <- function(x,
                       coords,
                       tree,
                       metric = c("mpd", "mntd"),
                       n_seeds = 50L,
                       method = "knn",
                       distance = c("euclidean", "haversine"),
                       parallel = TRUE,
                       n_cores = NULL,
                       progress = TRUE,
                       seed = NULL,
                       map = FALSE) {

  distance <- match.arg(distance)
  metric <- match.arg(metric, c("mpd", "mntd", "pd"), several.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    site_dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    if (progress) cli_info("Computing site distances")
    site_dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)

  # Get phylogenetic distance matrix
  if (inherits(tree, "phylo")) {
    check_suggests("ape")
    phylo_dist_mat <- as.matrix(ape::cophenetic.phylo(tree))

    # Match species order
    if (!is.null(colnames(x))) {
      if (!all(colnames(x) %in% tree$tip.label)) {
        stop("Some species in x are not in the tree")
      }
      phylo_dist_mat <- phylo_dist_mat[colnames(x), colnames(x)]
    }
  } else if (is.matrix(tree)) {
    phylo_dist_mat <- tree
    if ("pd" %in% metric) {
      warning("PD requires full tree object; removing 'pd' from metrics")
      metric <- setdiff(metric, "pd")
    }
  } else {
    stop("tree must be a phylo object or distance matrix")
  }

  stopifnot("Phylo distance matrix must match species" = ncol(phylo_dist_mat) == n_species)

  # Convert to presence/absence
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  # Remove pd for now (would need tree edge structure)
  cpp_metrics <- setdiff(metric, "pd")

  if (progress) cli_info(sprintf("Computing phylogenetic diversity (%s, %d seeds)",
                                  paste(metric, collapse = ", "), n_seeds))

  if (length(cpp_metrics) > 0) {
    result <- cpp_phylo_knn_parallel(species_pa, site_dist_mat, phylo_dist_mat,
                                      n_seeds, cpp_metrics, n_cores, progress)
  } else {
    result <- list()
  }

  # Handle PD separately if requested (needs R-side calculation)
  if ("pd" %in% metric && inherits(tree, "phylo")) {
    # Would need to implement PD accumulation
    # For now, skip or use picante
    warning("PD accumulation not yet implemented; use MPD/MNTD")
  }

  if (progress) cli_success("Done")

  # Compute per-site map values if requested
  site_values <- NULL
  if (map && length(cpp_metrics) > 0) {
    if (progress) cli_info("Computing per-site phylo map values (all sites as seeds)")
    map_result <- cpp_phylo_knn_parallel(species_pa, site_dist_mat, phylo_dist_mat,
                                          n_sites, cpp_metrics, n_cores, progress)

    site_values <- data.frame(
      site_id = seq_len(n_sites),
      x = coord_data$x,
      y = coord_data$y
    )
    for (i in seq_along(cpp_metrics)) {
      site_values[[cpp_metrics[i]]] <- map_result[[i]][cbind(seq_len(n_sites), n_sites)]
    }
    if (progress) cli_success("Map values computed")
  }

  structure(
    list(
      curves = result,
      metric = cpp_metrics,
      coords = coord_data,
      site_values = site_values,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_phylo"
  )
}


#' Spatial Functional Diversity Accumulation
#'
#' Compute spatial accumulation of functional diversity metrics based on traits.
#'
#' @param x A site-by-species matrix (abundance data recommended).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param traits A species-by-traits matrix. Row names should match species (columns of x).
#' @param metric Character vector. Metrics to compute:
#'   - `"fdis"`: Functional Dispersion (mean distance to centroid)
#'   - `"fric"`: Functional Richness (convex hull volume approximation)
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method. Default `"knn"`.
#' @param distance Character. Site distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final values for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#'
#' @return An object of class `spacc_func` containing:
#'   \item{curves}{Named list of matrices, one per metric (n_seeds x n_sites)}
#'   \item{metric}{Metrics computed}
#'   \item{coords, n_seeds, n_sites, method}{Parameters used}
#'
#' @details
#' Functional diversity metrics quantify trait space occupation:
#'
#' - **FDis (Functional Dispersion)**: Abundance-weighted mean distance from
#'   the community centroid in trait space. Captures functional divergence.
#'
#' - **FRic (Functional Richness)**: Volume of trait space occupied (convex hull).
#'   Requires more species than traits to compute.
#'
#' @references
#' Laliberté, E. & Legendre, P. (2010). A distance-based framework for measuring
#' functional diversity from multiple traits. Ecology, 91, 299-305.
#'
#' @seealso [FD::dbFD()] for comprehensive functional diversity analysis
#'
#' @examples
#' \dontrun{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 20, 2), nrow = 50)
#'
#' # Trait matrix (species x traits)
#' traits <- matrix(rnorm(20 * 3), nrow = 20)
#' rownames(traits) <- paste0("sp", 1:20)
#' colnames(species) <- rownames(traits)
#'
#' func <- spaccFunc(species, coords, traits, metric = c("fdis", "fric"))
#' plot(func)
#' }
#'
#' @export
spaccFunc <- function(x,
                      coords,
                      traits,
                      metric = c("fdis", "fric"),
                      n_seeds = 50L,
                      method = "knn",
                      distance = c("euclidean", "haversine"),
                      parallel = TRUE,
                      n_cores = NULL,
                      progress = TRUE,
                      seed = NULL,
                      map = FALSE) {

  distance <- match.arg(distance)
  metric <- match.arg(metric, c("fdis", "fric"), several.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  traits <- as.matrix(traits)

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    site_dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    if (progress) cli_info("Computing site distances")
    site_dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)
  n_traits <- ncol(traits)

  # Match traits to species
  if (!is.null(colnames(x)) && !is.null(rownames(traits))) {
    if (!all(colnames(x) %in% rownames(traits))) {
      stop("Some species in x are not in traits")
    }
    traits <- traits[colnames(x), , drop = FALSE]
  }

  stopifnot("Traits must have one row per species" = nrow(traits) == n_species)

  # Keep abundance data
  storage.mode(x) <- "integer"

  if (progress) cli_info(sprintf("Computing functional diversity (%s, %d seeds)",
                                  paste(metric, collapse = ", "), n_seeds))

  result <- cpp_func_knn_parallel(x, site_dist_mat, traits,
                                   n_seeds, metric, n_cores, progress)

  if (progress) cli_success("Done")

  # Compute per-site map values if requested
  site_values <- NULL
  if (map) {
    if (progress) cli_info("Computing per-site functional map values (all sites as seeds)")
    map_result <- cpp_func_knn_parallel(x, site_dist_mat, traits,
                                         n_sites, metric, n_cores, progress)

    site_values <- data.frame(
      site_id = seq_len(n_sites),
      x = coord_data$x,
      y = coord_data$y
    )
    for (i in seq_along(metric)) {
      site_values[[metric[i]]] <- map_result[[i]][cbind(seq_len(n_sites), n_sites)]
    }
    if (progress) cli_success("Map values computed")
  }

  structure(
    list(
      curves = result,
      metric = metric,
      traits = traits,
      coords = coord_data,
      site_values = site_values,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      n_traits = n_traits,
      method = method,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_func"
  )
}


# S3 methods for spacc_phylo -----------------------------------------------

#' @export
print.spacc_phylo <- function(x, ...) {
  cat(sprintf("spacc phylogenetic diversity: %d sites, %d species, %d seeds\n",
              x$n_sites, x$n_species, x$n_seeds))
  cat(sprintf("Metrics: %s\n", paste(x$metric, collapse = ", ")))
  invisible(x)
}


#' @export
summary.spacc_phylo <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2

  result <- lapply(seq_along(object$metric), function(i) {
    metric_name <- object$metric[i]
    mat <- object$curves[[i]]

    data.frame(
      metric = metric_name,
      sites = 1:object$n_sites,
      mean = colMeans(mat),
      lower = apply(mat, 2, stats::quantile, alpha),
      upper = apply(mat, 2, stats::quantile, 1 - alpha)
    )
  })

  do.call(rbind, result)
}


#' @export
plot.spacc_phylo <- function(x, type = c("curve", "map"), ci = TRUE, ci_alpha = 0.2,
                              metric = NULL, point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$site_values)) stop("No map data. Rerun spaccPhylo() with map = TRUE.")
    if (is.null(metric)) metric <- x$metric[1]
    stopifnot("metric not found" = metric %in% names(x$site_values))
    return(plot_spatial_map(x$site_values, metric,
                            title = sprintf("Phylogenetic diversity map: %s", toupper(metric)),
                            subtitle = sprintf("%d sites, %s method", x$n_sites, x$method),
                            point_size = point_size, palette = palette))
  }

  check_suggests("ggplot2")

  summ <- summary(x)
  summ$metric <- factor(toupper(summ$metric))

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = sites, y = mean, color = metric))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = metric),
      alpha = ci_alpha, color = NA
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "Sites accumulated",
      y = "Phylogenetic diversity",
      color = "Metric",
      fill = "Metric",
      title = "Spatial Phylogenetic Diversity Accumulation"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}


#' @export
as_sf.spacc_phylo <- function(x, crs = NULL) {
  if (is.null(x$site_values)) stop("No map data. Rerun spaccPhylo() with map = TRUE.")
  as_sf_from_df(x$site_values, crs = crs)
}


# S3 methods for spacc_func -----------------------------------------------

#' @export
print.spacc_func <- function(x, ...) {
  cat(sprintf("spacc functional diversity: %d sites, %d species, %d traits, %d seeds\n",
              x$n_sites, x$n_species, x$n_traits, x$n_seeds))
  cat(sprintf("Metrics: %s\n", paste(x$metric, collapse = ", ")))
  invisible(x)
}


#' @export
summary.spacc_func <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2

  result <- lapply(seq_along(object$metric), function(i) {
    metric_name <- object$metric[i]
    mat <- object$curves[[i]]

    data.frame(
      metric = metric_name,
      sites = 1:object$n_sites,
      mean = colMeans(mat),
      lower = apply(mat, 2, stats::quantile, alpha),
      upper = apply(mat, 2, stats::quantile, 1 - alpha)
    )
  })

  do.call(rbind, result)
}


#' @export
plot.spacc_func <- function(x, type = c("curve", "map"), ci = TRUE, ci_alpha = 0.2,
                             metric = NULL, point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$site_values)) stop("No map data. Rerun spaccFunc() with map = TRUE.")
    if (is.null(metric)) metric <- x$metric[1]
    stopifnot("metric not found" = metric %in% names(x$site_values))
    return(plot_spatial_map(x$site_values, metric,
                            title = sprintf("Functional diversity map: %s", toupper(metric)),
                            subtitle = sprintf("%d sites, %s method", x$n_sites, x$method),
                            point_size = point_size, palette = palette))
  }

  check_suggests("ggplot2")

  summ <- summary(x)
  summ$metric <- factor(toupper(summ$metric))

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = sites, y = mean, color = metric))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = metric),
      alpha = ci_alpha, color = NA
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "Sites accumulated",
      y = "Functional diversity",
      color = "Metric",
      fill = "Metric",
      title = "Spatial Functional Diversity Accumulation"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}


#' @export
as_sf.spacc_func <- function(x, crs = NULL) {
  if (is.null(x$site_values)) stop("No map data. Rerun spaccFunc() with map = TRUE.")
  as_sf_from_df(x$site_values, crs = crs)
}
