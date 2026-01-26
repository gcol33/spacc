#' Spatial Species Accumulation Curves
#'
#' Compute species accumulation curves using spatially-explicit nearest-neighbor
#' sampling with C++ backend for performance.
#'
#' @param x A site-by-species matrix (rows = sites, cols = species) with
#'   presence/absence (0/1) or abundance data. Can also be a data.frame.
#' @param coords A data.frame with columns `x` and `y` containing site coordinates,
#'   or a `spacc_dist` object from [distances()].
#' @param n_seeds Integer. Number of random starting points for uncertainty
#'   quantification. Default 50.
#' @param method Character. One of `"knn"` (k-Nearest Neighbor), `"kncn"`
#'   (k-Nearest Centroid Neighbor), or `"random"` (random order, null model).
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`
#'   (for lat/lon coordinates). Ignored if `coords` is a `spacc_dist` object.
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
#' # With options
#' sac <- spacc(species, coords,
#'              n_seeds = 100,
#'              method = "knn",
#'              distance = "haversine",
#'              progress = TRUE)
#'
#' # Reuse distance matrix
#' d <- distances(coords, method = "haversine")
#' sac1 <- spacc(native_species, d)
#' sac2 <- spacc(alien_species, d)
#' }
#'
#' @export
spacc <- function(x,
                  coords,
                  n_seeds = 50L,
                  method = c("knn", "kncn", "random"),
                  distance = c("euclidean", "haversine"),
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
    dist_mat <- coords
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

  # Compute distance matrix if not provided
  if (is.null(dist_mat)) {
    if (progress) cli_info(sprintf("Computing distances (%d x %d)", n_sites, n_sites))
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  # Run accumulation curves
  if (progress) cli_info(sprintf("Running %s accumulation (%d seeds, %d cores)", method, n_seeds, n_cores))

  curves <- switch(method,
    knn = cpp_knn_parallel(species_pa, dist_mat, n_seeds, n_cores, progress),
    kncn = cpp_kncn_parallel(species_pa, coord_data$x, coord_data$y, n_seeds, n_cores, progress),
    random = cpp_random_parallel(species_pa, n_seeds, n_cores, progress)
  )

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
      call = match.call()
    ),
    class = "spacc"
  )
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
