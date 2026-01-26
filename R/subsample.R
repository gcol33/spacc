#' Spatial Subsampling
#'
#' Reduce spatial autocorrelation by subsampling sites using various methods.
#'
#' @param coords A data.frame with columns `x` and `y` containing site coordinates.
#' @param n Integer. Target number of sites to retain. If `NULL`, determined
#'   by `cell_size` or `min_dist`.
#' @param method Character. Subsampling method: `"grid"` (default), `"random"`,
#'   or `"thinning"`.
#' @param cell_size Numeric. Grid cell size for `method = "grid"`. One site
#'   retained per cell.
#' @param min_dist Numeric. Minimum distance between retained sites for
#'   `method = "thinning"`.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return Integer vector of row indices to retain.
#'
#' @details
#' **Methods:**
#' - `"grid"`: Overlay a grid and retain one random site per cell.
#' - `"random"`: Simple random subsample of n sites.
#' - `"thinning"`: Iteratively remove sites until minimum distance is achieved.
#'
#' @examples
#' \dontrun{
#' # Grid-based subsampling
#' keep <- subsample(coords, method = "grid", cell_size = 10)
#' sac <- spacc(species[keep, ], coords[keep, ])
#'
#' # Minimum distance thinning
#' keep <- subsample(coords, method = "thinning", min_dist = 5)
#' }
#'
#' @export
subsample <- function(coords,
                      n = NULL,
                      method = c("grid", "random", "thinning"),
                      cell_size = NULL,
                      min_dist = NULL,
                      seed = NULL) {

  method <- match.arg(method)

  stopifnot(
    "coords must have x and y columns" = all(c("x", "y") %in% names(coords))
  )

  if (!is.null(seed)) set.seed(seed)

  n_sites <- nrow(coords)
  x <- coords$x
  y <- coords$y

  if (method == "random") {
    stopifnot("n must be specified for random subsampling" = !is.null(n))
    return(sample(n_sites, min(n, n_sites)))
  }

  if (method == "grid") {
    if (is.null(cell_size)) {
      stopifnot("Either cell_size or n must be specified" = !is.null(n))
      # Estimate cell size to get approximately n sites
      extent_x <- diff(range(x))
      extent_y <- diff(range(y))
      area <- extent_x * extent_y
      cell_size <- sqrt(area / n)
    }

    # Assign each site to a grid cell
    cell_x <- floor(x / cell_size)
    cell_y <- floor(y / cell_size)
    cell_id <- paste(cell_x, cell_y, sep = "_")

    # Sample one site per cell
    cells <- unique(cell_id)
    keep <- vapply(cells, function(cell) {
      idx <- which(cell_id == cell)
      sample(idx, 1)
    }, integer(1))

    return(unname(keep))
  }

  if (method == "thinning") {
    stopifnot("min_dist must be specified for thinning" = !is.null(min_dist))

    # Greedy thinning: iteratively remove closest pairs
    keep <- seq_len(n_sites)
    dist_mat <- as.matrix(stats::dist(cbind(x, y)))

    repeat {
      if (length(keep) <= 1) break

      sub_dist <- dist_mat[keep, keep, drop = FALSE]
      diag(sub_dist) <- Inf

      min_d <- min(sub_dist)
      if (min_d >= min_dist) break

      # Find the pair with minimum distance
      idx <- which(sub_dist == min_d, arr.ind = TRUE)[1, ]

      # Remove the one with more close neighbors
      neighbors1 <- sum(sub_dist[idx[1], ] < min_dist)
      neighbors2 <- sum(sub_dist[idx[2], ] < min_dist)

      remove_idx <- ifelse(neighbors1 >= neighbors2, idx[1], idx[2])
      keep <- keep[-remove_idx]
    }

    return(keep)
  }
}
