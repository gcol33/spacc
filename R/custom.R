#' Spatial Accumulation of a Custom Diversity Metric
#'
#' Accumulate any user-supplied diversity index along a spatial ordering of
#' sites. At each accumulation step the cumulative community is passed to
#' `fun`, which returns a single number. This is the general escape hatch
#' behind the built-in metric functions: use it for indices that spacc does
#' not implement directly.
#'
#' @param x A site-by-species matrix (rows = sites, cols = species), abundance
#'   or presence/absence.
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param fun A function applied to the cumulative community at each step. It
#'   receives a named numeric vector of length `ncol(x)` (cumulative summed
#'   abundances, or 0/1 incidences when `incidence = TRUE`) plus any arguments
#'   passed through `...`, and must return a single numeric value.
#' @param ... Additional arguments passed to `fun`.
#' @param method Character. Spatial ordering of sites: `"knn"` (default),
#'   `"kncn"`, `"random"`, `"radius"`, or `"collector"`.
#' @param incidence Logical. If `TRUE`, `fun` receives 0/1 incidences instead
#'   of summed abundances. Default `FALSE`.
#' @param n_seeds Integer. Number of random starting points / orderings.
#'   Ignored for `"collector"` (a single data-order curve). Default 50.
#' @param distance Character. `"euclidean"` or `"haversine"`.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return An object of class `spacc_diversity` that inherits from `spacc`, so
#'   the standard `summary()`, `plot()`, `as.data.frame()` and `predict()`
#'   methods apply. `curves` is an `n_seeds x n_sites` matrix of the metric
#'   along the accumulation.
#'
#' @details
#' The site ordering reuses the same spatial traversals as the built-in
#' methods (nearest-neighbour, nearest-centroid, random, distance-rank, or
#' data order), then evaluates `fun` on the accumulating community. Because the
#' index is an arbitrary R function, this trades the speed of the compiled
#' metrics for full flexibility.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(40), y = runif(40))
#' species <- matrix(rpois(40 * 20, 2), nrow = 40)
#'
#' # Shannon entropy along the accumulation
#' shannon <- function(comm) {
#'   p <- comm[comm > 0] / sum(comm)
#'   -sum(p * log(p))
#' }
#' div <- spaccDiversity(species, coords, shannon, n_seeds = 20)
#' plot(div)
#' }
#'
#' @seealso [spaccHill()], [spaccPhylo()], [spaccFunc()] for built-in metrics.
#'
#' @export
spaccDiversity <- function(x, coords, fun, ...,
                           method = c("knn", "kncn", "random", "radius", "collector"),
                           incidence = FALSE,
                           n_seeds = 50L,
                           distance = c("euclidean", "haversine"),
                           progress = TRUE,
                           seed = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  stopifnot("fun must be a function" = is.function(fun))

  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)

  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    dist_mat <- NULL
  }

  stopifnot(
    "x and coords must have same number of rows" = nrow(x) == nrow(coord_data),
    "n_seeds must be positive" = n_seeds > 0
  )

  n_sites <- nrow(x)
  n_species <- ncol(x)
  sp_names <- colnames(x)

  comm_mat <- if (incidence) (x > 0) * 1 else x
  storage.mode(comm_mat) <- "double"

  if (method %in% c("knn", "radius") && is.null(dist_mat)) {
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  # Build accumulation orders (1-based), one row per ordering
  if (method == "collector") {
    orders <- matrix(seq_len(n_sites), nrow = 1L)
  } else {
    seeds0 <- as.integer(sample(n_sites, n_seeds, replace = TRUE) - 1L)
    orders <- switch(method,
      knn    = cpp_knn_order(dist_mat, seeds0) + 1L,
      kncn   = cpp_kncn_order(coord_data$x, coord_data$y, seeds0) + 1L,
      random = t(vapply(seq_len(n_seeds), function(i) sample.int(n_sites), integer(n_sites))),
      radius = t(vapply(seeds0 + 1L, function(s) order(dist_mat[s, ]), integer(n_sites)))
    )
  }
  n_curves <- nrow(orders)

  if (progress) cli_info(sprintf("Running custom diversity accumulation (%s, %d ordering%s)",
                                  method, n_curves, if (n_curves == 1L) "" else "s"))

  curves <- matrix(NA_real_, nrow = n_curves, ncol = n_sites)
  for (s in seq_len(n_curves)) {
    ord <- orders[s, ]
    cum <- numeric(n_species)
    for (k in seq_len(n_sites)) {
      cum <- cum + comm_mat[ord[k], ]
      val <- fun(stats::setNames(cum, sp_names), ...)
      if (length(val) != 1L || !is.numeric(val)) {
        stop("`fun` must return a single numeric value", call. = FALSE)
      }
      curves[s, k] <- val
    }
  }

  if (progress) cli_success("Done")

  structure(
    list(
      curves = curves,
      coords = coord_data,
      n_seeds = n_curves,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      distance = distance,
      fun = fun,
      call = match.call()
    ),
    class = c("spacc_diversity", "spacc")
  )
}


#' @export
print.spacc_diversity <- function(x, ...) {
  cat(sprintf("spacc custom diversity: %d sites, %d species, %d ordering%s (%s)\n",
              x$n_sites, x$n_species, x$n_seeds,
              if (x$n_seeds == 1L) "" else "s", x$method))
  invisible(x)
}


#' @export
plot.spacc_diversity <- function(x, ...,
                                 ylab = "Cumulative diversity",
                                 title = "Custom Diversity Accumulation") {
  check_suggests("ggplot2")
  # Reuse the spacc accumulation plot, relabelling for an arbitrary metric.
  plot.spacc(x, ...) + ggplot2::labs(y = ylab, title = title)
}
