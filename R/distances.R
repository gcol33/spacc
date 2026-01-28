#' Compute Distance Matrix
#'
#' Pre-compute pairwise distances between sites for reuse across multiple
#' `spacc()` calls. Supports sf objects with accurate geodesic distances
#' for global-scale studies.
#'
#' @param x Site locations. Can be:
#'   - A data.frame with columns `x` and `y`
#'   - An sf object with POINT geometries
#'   - An sfc_POINT object
#' @param method Character. Distance method:
#'   - `"euclidean"`: Euclidean distance (for projected coordinates)
#'   - `"haversine"`: Great-circle distance (for lat/lon, fast approximation
#'   - `"geodesic"`: Accurate ellipsoidal distance via sf/S2 (for global scale)
#'   Default is auto-detected from CRS when using sf objects.
#' @param fun Optional custom distance function. Must take two coordinate
#'   vectors (x, y) and return a distance matrix. Overrides `method`.
#' @param which For sf objects, column name containing the geometry.
#'   Default uses active geometry.
#'
#' @return An object of class `spacc_dist` containing the distance matrix
#'   with coordinates stored as an attribute.
#'
#' @details
#' For continental and global-scale studies, use sf objects with geographic
#' CRS (e.g., EPSG:4326). The function will automatically use accurate
#' geodesic distances via the S2 spherical geometry library.
#'
#' For smaller study areas with projected coordinates (UTM, etc.), Euclidean
#' distance is appropriate and faster.
#'
#' @examples
#' \dontrun{
#' # Data frame with coordinates
#' d <- distances(coords, method = "haversine")
#'
#' # sf object - auto-detects appropriate method
#' library(sf)
#' pts_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
#' d <- distances(pts_sf)  # uses geodesic automatically
#'
#' # Force geodesic for maximum accuracy
#' d <- distances(pts_sf, method = "geodesic")
#'
#' # Reuse for multiple analyses
#' sac_native <- spacc(native_species, d)
#' sac_alien <- spacc(alien_species, d)
#' }
#'
#' @export
distances <- function(x, method = NULL, fun = NULL, which = NULL) {
  UseMethod("distances")
}


#' @export
distances.data.frame <- function(x, method = NULL, fun = NULL, which = NULL) {

  if (is.null(method)) method <- "euclidean"
  method <- match.arg(method, c("euclidean", "haversine"))

  stopifnot(
    "coords must have x and y columns" = all(c("x", "y") %in% names(x))
  )

  if (!is.null(fun)) {
    dist_mat <- fun(x$x, x$y)
  } else {
    dist_mat <- cpp_distance_matrix(
      as.numeric(x$x),
      as.numeric(x$y),
      method
    )
  }

  structure(
    dist_mat,
    coords = x[, c("x", "y"), drop = FALSE],
    method = method,
    crs = NA,
    class = c("spacc_dist", "matrix", "array")
  )
}


#' @export
distances.sf <- function(x, method = NULL, fun = NULL, which = NULL) {

  check_suggests("sf")

  # Extract geometry
  if (!is.null(which)) {
    geom <- x[[which]]
  } else {
    geom <- sf::st_geometry(x)
  }

  # Check geometry type
  geom_type <- unique(sf::st_geometry_type(geom))
  if (!all(geom_type %in% c("POINT", "MULTIPOINT"))) {
    stop("Only POINT geometries are supported. Got: ", paste(geom_type, collapse = ", "))
  }

  # Get CRS info
 crs <- sf::st_crs(x)
  is_geographic <- !is.na(crs) && sf::st_is_longlat(x)

  # Auto-detect method if not specified
  if (is.null(method)) {
    method <- if (is_geographic) "geodesic" else "euclidean"
    message(sprintf("Auto-detected CRS: using '%s' distance", method))
  }
  method <- match.arg(method, c("euclidean", "haversine", "geodesic"))

  # Extract coordinates
  coords_mat <- sf::st_coordinates(geom)
  coord_df <- data.frame(x = coords_mat[, 1], y = coords_mat[, 2])

  if (!is.null(fun)) {
    dist_mat <- fun(coord_df$x, coord_df$y)
  } else if (method == "geodesic") {
    # Use sf::st_distance for accurate geodesic distances
    dist_obj <- sf::st_distance(geom)
    dist_mat <- matrix(as.numeric(dist_obj), nrow = nrow(coords_mat))
  } else {
    # Use fast C++ implementation
    dist_mat <- cpp_distance_matrix(
      as.numeric(coord_df$x),
      as.numeric(coord_df$y),
      method
    )
  }

  structure(
    dist_mat,
    coords = coord_df,
    sf_geometry = geom,
    method = method,
    crs = crs,
    class = c("spacc_dist", "matrix", "array")
  )
}


#' @export
distances.sfc <- function(x, method = NULL, fun = NULL, which = NULL) {
  # Convert sfc to sf and dispatch
  sf_obj <- sf::st_as_sf(data.frame(id = seq_along(x)), geometry = x)
  distances.sf(sf_obj, method = method, fun = fun, which = NULL)
}


#' @export
print.spacc_dist <- function(x, ...) {
  n <- nrow(x)
  method <- attr(x, "method")
  crs <- attr(x, "crs")

  cat(sprintf("spacc distance matrix: %d x %d\n", n, n))
  cat(sprintf("Method: %s\n", method))

  if (!is.null(crs) && is.list(crs)) {
    if (!is.na(crs$epsg)) {
      cat(sprintf("CRS: EPSG:%s\n", crs$epsg))
    } else if (!is.na(crs$input)) {
      cat(sprintf("CRS: %s\n", substr(crs$input, 1, 50)))
    }
  }

  # Distance summary
  d_vals <- x[lower.tri(x)]
  cat(sprintf("Distance range: %.1f - %.1f", min(d_vals), max(d_vals)))
  if (method == "geodesic") cat(" (meters)")
  cat("\n")

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
        x = sprintf("Distance (%s)", attr(x, "method")),
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
