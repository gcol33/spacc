#' Compute Distance Matrix
#'
#' Pre-compute pairwise distances between sites for reuse across multiple
#' `spacc()` calls. Useful when comparing groups (e.g., native vs alien)
#' on the same set of sites.
#'
#' @param coords A data.frame with columns `x` and `y` containing site coordinates.
#' @param method Character. Distance method: `"euclidean"` (default) or
#'   `"haversine"` (for lat/lon coordinates in degrees).
#' @param fun Optional custom distance function. Must take two coordinate
#'   vectors and return a distance matrix. Overrides `method`.
#'
#' @return An object of class `spacc_dist` containing the distance matrix
#'   with coordinates stored as an attribute.
#'
#' @examples
#' \dontrun{
#' # Pre-compute distances
#' d <- distances(coords, method = "haversine")
#'
#' # Reuse for multiple analyses
#' sac_native <- spacc(native_species, d)
#' sac_alien <- spacc(alien_species, d)
#' }
#'
#' @export
distances <- function(coords, method = c("euclidean", "haversine"), fun = NULL) {

  method <- match.arg(method)

  stopifnot(
    "coords must have x and y columns" = all(c("x", "y") %in% names(coords))
  )

  if (!is.null(fun)) {
    dist_mat <- fun(coords$x, coords$y)
  } else {
    dist_mat <- cpp_distance_matrix(
      as.numeric(coords$x),
      as.numeric(coords$y),
      method
    )
  }

  structure(
    dist_mat,
    coords = coords,
    method = method,
    class = c("spacc_dist", "matrix", "array")
  )
}


#' @export
print.spacc_dist <- function(x, ...) {
  n <- nrow(x)
  method <- attr(x, "method")
  cat(sprintf("spacc distance matrix: %d x %d (%s)\n", n, n, method))
  invisible(x)
}


#' @export
plot.spacc_dist <- function(x, type = c("histogram", "heatmap"), ...) {
  check_suggests("ggplot2")

  type <- match.arg(type)

  d_vals <- x[lower.tri(x)]

  if (type == "histogram") {
    df <- data.frame(distance = d_vals)
    ggplot2::ggplot(df, ggplot2::aes(x = distance)) +
      ggplot2::geom_histogram(fill = "#4CAF50", color = "white", bins = 50) +
      ggplot2::labs(
        title = "Pairwise Distance Distribution",
        x = "Distance",
        y = "Count"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  } else {
    df <- expand.grid(i = seq_len(nrow(x)), j = seq_len(ncol(x)))
    df$distance <- as.vector(x)
    ggplot2::ggplot(df, ggplot2::aes(x = i, y = j, fill = distance)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(
        title = "Distance Matrix",
        x = "Site",
        y = "Site"
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 12)
  }
}
