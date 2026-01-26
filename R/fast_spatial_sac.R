#' Fast Spatial Species Accumulation Curve
#'
#' Compute species accumulation curves using spatially-explicit nearest-neighbor
#' sampling with C++ backend for performance.
#'
#' @param species A site-by-species matrix (rows = sites, cols = species) with
#'   presence/absence (0/1) or abundance data. Can also be a data.frame.
#' @param coords A data.frame with columns `x` and `y` containing site coordinates.
#'   Must have same number of rows as `species`.
#' @param n_seeds Integer. Number of random starting points for uncertainty
#'   quantification. Default 50.
#' @param method Character. Either `"knn"` (k-Nearest Neighbor) or `"kncn"`
#'
#' @return An object of class `"spatial_sac"` containing:
#'   \item{curves}{Matrix of cumulative species counts (n_seeds x n_sites)}
#'   \item{coords}{Original coordinates}
#'   \item{n_seeds}{Number of seeds used
#'   \item{method}{Method used}
#'   \item{n_species}{Total species in dataset}
#'
#' @examples
#' \dontrun
#' # Create example data
#' set.seed(42)
#' n_sites <- 100
#' n_species <- 50
#'
#' coords <- data.frame(
#'   x = runif(n_sites, 0, 100),
#'   y = runif(n_sites, 0, 100
#' )
#'
#' species <- matrix(
#'   rbinom(n_sites * n_species, 1, 0.3),
#'   nrow = n_sites
#' )
#'
#' sac <- fast_spatial_sac(species, coords, n_seeds = 20)
#' plot(sac)
#' }
#'
#' @export
fast_spatial_sac <- function(species, coords, n_seeds = 50L, method = c("knn", "kncn")) {

 method <- match.arg(method)

 # Input validation
 species <- as.matrix(species)
 stopifnot(
   "coords must have x and y columns" = all(c("x", "y") %in% names(coords)),
   "species and coords must have same number of rows" = nrow(species) == nrow(coords),
   "n_seeds must be positive integer" = n_seeds > 0
 )

 n_sites <- nrow(species)
 n_species_total <- ncol(species)

 # Convert to presence/absence if abundance
 species_pa <- (species > 0) * 1L

 # Compute distance matrix (C++)
 dist_mat <- fast_distance_matrix(coords$x, coords$y)

 # Run accumulation curves (C++ parallel)
 if (method == "knn") {
   curves <- cpp_knn_parallel(species_pa, dist_mat, n_seeds)
 } else {
   curves <- cpp_kncn_parallel(species_pa, dist_mat, n_seeds)
 }

 structure(
   list(
     curves = curves,
     coords = coords,
     n_seeds = n_seeds,
     method = method,
     n_sites = n_sites,
     n_species = n_species_total,
     call = match.call()
   ),
   class = "spatial_sac"
 )
}


#' Fast Distance Matrix Computation
#'
#' Compute pairwise Euclidean distances using C++.
#'
#' @param x Numeric vector of x coordinates
#' @param y Numeric vector of y coordinates
#'
#' @return Symmetric distance matrix (n x n)
#'
#' @export
fast_distance_matrix <- function(x, y) {
 stopifnot(length(x) == length(y))
 cpp_distance_matrix(as.numeric(x), as.numeric(y))
}
