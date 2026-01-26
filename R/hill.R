#' Spatial Accumulation with Hill Numbers
#'
#' Compute spatial species accumulation curves using Hill numbers (effective
#' number of species) instead of raw richness. Hill numbers unify diversity
#' measures: q=0 is richness, q=1 is exponential Shannon, q=2 is inverse Simpson.
#'
#' @param x A site-by-species matrix (rows = sites, cols = species) with
#'   presence/absence (0/1) or abundance data.
#' @param coords A data.frame with columns `x` and `y` containing site coordinates,
#'   or a `spacc_dist` object from [distances()].
#' @param q Numeric vector. Orders of diversity to compute. Default `c(0, 1, 2)`.
#'   - q = 0: Species richness
#'   - q = 1: Exponential of Shannon entropy (effective common species)
#'   - q = 2: Inverse Simpson (effective dominant species)
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method: `"knn"` (default).
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores. Default `NULL` uses `detectCores() - 1`.
#' @param progress Logical. Show progress bar? Default `TRUE`.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return An object of class `spacc_hill` containing:
#'   \item{curves}{Named list of matrices, one per q value (n_seeds x n_sites)}
#'   \item{q}{Vector of q values used}
#'   \item{coords}{Original coordinates}
#'   \item{n_seeds}{Number of seeds}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Total species}
#'   \item{method}{Method used}
#'
#' @details
#' Hill numbers (Chao et al. 2014) provide a unified framework for diversity
#' measurement. Unlike raw richness (q=0), higher-order Hill numbers (q=1, q=2)
#' down-weight rare species, providing different perspectives on diversity.
#'
#' The spatial accumulation of Hill numbers can reveal scale-dependent diversity
#' patterns missed by richness alone.
#'
#' @references
#' Chao, A., Gotelli, N.J., Hsieh, T.C., Sander, E.L., Ma, K.H., Colwell, R.K.
#' & Ellison, A.M. (2014). Rarefaction and extrapolation with Hill numbers:
#' a framework for sampling and estimation in species diversity studies.
#' Ecological Monographs, 84, 45-67.
#'
#' @seealso [spacc()] for richness-only accumulation, [iNEXT::iNEXT()] for
#'   non-spatial Hill number rarefaction
#'
#' @examples
#' \dontrun{
#' # Compare diversity at different orders
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' hill <- spaccHill(species, coords, q = c(0, 1, 2))
#' plot(hill)
#'
#' # Extract summary at final site
#' summary(hill)
#' }
#'
#' @export
spaccHill <- function(x,
                      coords,
                      q = c(0, 1, 2),
                      n_seeds = 50L,
                      method = "knn",
                      distance = c("euclidean", "haversine"),
                      parallel = TRUE,
                      n_cores = NULL,
                      progress = TRUE,
                      seed = NULL) {

  distance <- match.arg(distance)

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  # Validate inputs
  x <- as.matrix(x)
  stopifnot(
    "q must be numeric and >= 0" = is.numeric(q) && all(q >= 0),
    "n_seeds must be positive" = n_seeds > 0
  )

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot(
      "coords must have x and y columns" = all(c("x", "y") %in% names(coords))
    )
    coord_data <- coords
    if (progress) cli_info(sprintf("Computing distances (%d x %d)", nrow(x), nrow(x)))
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  stopifnot(
    "x and coords must have same rows" = nrow(x) == nrow(coord_data)
  )

  n_sites <- nrow(x)
  n_species <- ncol(x)

  # Ensure integer matrix (abundances preserved)
  storage.mode(x) <- "integer"

  if (progress) cli_info(sprintf("Computing Hill numbers (q = %s, %d seeds)",
                                  paste(q, collapse = ", "), n_seeds))

  # Call C++ function
  curves <- cpp_knn_hill_parallel(x, dist_mat, n_seeds, q, n_cores, progress)

  if (progress) cli_success("Done")

  structure(
    list(
      curves = curves,
      q = q,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_hill"
  )
}


#' @export
print.spacc_hill <- function(x, ...) {
  cat(sprintf("spacc Hill numbers: %d sites, %d species, %d seeds\n",
              x$n_sites, x$n_species, x$n_seeds))
  cat(sprintf("Orders (q): %s\n", paste(x$q, collapse = ", ")))
  cat(sprintf("Method: %s\n", x$method))
  invisible(x)
}


#' @export
summary.spacc_hill <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2

  result <- lapply(seq_along(object$q), function(i) {
    q_val <- object$q[i]
    q_name <- names(object$curves)[i]
    mat <- object$curves[[i]]

    mean_curve <- colMeans(mat)
    lower <- apply(mat, 2, stats::quantile, probs = alpha)
    upper <- apply(mat, 2, stats::quantile, probs = 1 - alpha)

    data.frame(
      q = q_val,
      sites = 1:object$n_sites,
      mean = mean_curve,
      lower = lower,
      upper = upper
    )
  })

  do.call(rbind, result)
}


#' @export
plot.spacc_hill <- function(x, ci = TRUE, ci_alpha = 0.2, ...) {
  check_suggests("ggplot2")

  summ <- summary(x)
  summ$q_label <- factor(paste0("q = ", summ$q),
                         levels = paste0("q = ", sort(unique(summ$q))))

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = sites, y = mean, color = q_label))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = q_label),
      alpha = ci_alpha, color = NA
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "Sites accumulated",
      y = "Hill number (effective species)",
      color = "Order",
      fill = "Order",
      title = "Spatial Hill Number Accumulation"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
