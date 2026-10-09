#' Spatial Hill Numbers at Standardized Coverage
#'
#' Compute spatial accumulation of Hill numbers alongside sample coverage,
#' enabling standardized comparison at equal completeness levels (iNEXT-style
#' analysis applied spatially).
#'
#' @param x A site-by-species matrix (rows = sites, cols = species) with
#'   abundance data.
#' @param coords A data.frame with columns `x` and `y` containing site
#'   coordinates, or a `spacc_dist` object from [distances()].
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#' @param target_coverage Numeric vector. Target coverage levels for
#'   standardization. Default `NULL` (no standardization). Values in (0, 1).
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores. Default `NULL` uses all minus one.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed for reproducibility.
#' @param method Character. Accumulation method: `"knn"` or `"nn_walk"`.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#'
#' @return An object of class `spacc_hill_coverage` containing:
#'   \item{hills}{Named list of n_seeds x n_sites matrices (one per q)}
#'   \item{coverage}{n_seeds x n_sites matrix of Chao-Jost coverage}
#'   \item{pooled}{Species abundances pooled over all sites}
#'   \item{standardized}{List of per-q values interpolated at target coverage,
#'     or `NULL` if `target_coverage` is `NULL`}
#'   \item{q}{Vector of q values}
#'   \item{target_coverage}{Target coverage levels used}
#'   \item{coords}{Original coordinates}
#'   \item{n_seeds}{Number of seeds}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Total species}
#'
#' @details
#' Combines spatial kNN accumulation with simultaneous Hill number and
#' coverage computation. At each accumulation step, both the Chao-Jost
#' sample coverage and Hill numbers for all requested orders are calculated.
#'
#' When `target_coverage` is specified, Hill numbers are interpolated at
#' those coverage levels using the existing `interpolate_at_coverage()` C++
#' function, enabling fair comparison between sites with different sampling
#' completeness.
#'
#' @references
#' Chao, A. & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. Ecology, 93,
#' 2533-2547.
#'
#' Chao, A., Gotelli, N.J., Hsieh, T.C., Sander, E.L., Ma, K.H., Colwell,
#' R.K. & Ellison, A.M. (2014). Rarefaction and extrapolation with Hill
#' numbers. Ecological Monographs, 84, 45-67.
#'
#' @seealso [spaccHill()] for Hill accumulation without coverage,
#'   [spaccCoverage()] for coverage accumulation without Hill numbers
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(40), y = runif(40))
#' species <- matrix(rpois(40 * 20, 2), nrow = 40)
#'
#' hc <- spaccHillCoverage(species, coords, n_seeds = 10, progress = FALSE)
#' plot(hc)
#'
#' # Standardize at 90% coverage
#' hc2 <- spaccHillCoverage(species, coords, target_coverage = 0.9,
#'                           n_seeds = 10, progress = FALSE)
#' summary(hc2)
#' }
#'
#' @export
spaccHillCoverage <- function(x,
                               coords,
                               q = c(0, 1, 2),
                               target_coverage = NULL,
                               n_seeds = 50L,
                               distance = c("euclidean", "haversine"),
                               parallel = TRUE,
                               n_cores = NULL,
                               progress = TRUE,
                               seed = NULL,
                               method = c("knn", "nn_walk"),
                               focal_points = NULL,
                               focal_domain = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  if (!is.null(seed)) set.seed(seed)
  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  stopifnot(
    "q must be numeric and >= 0" = is.numeric(q) && all(q >= 0),
    "n_seeds must be positive" = n_seeds > 0
  )

  if (!is.null(target_coverage)) {
    stopifnot(
      "target_coverage must be in (0, 1)" = all(target_coverage > 0 & target_coverage < 1)
    )
  }

  site <- .resolve_coords(coords)
  coord_data <- site$coords
  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)
  storage.mode(x) <- "integer"

  if (progress) cli_info(sprintf("Computing Hill numbers + coverage (q = %s, %d seeds)",
                                  paste(q, collapse = ", "), n_seeds))

  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   site$dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)
  n_seeds <- ordering$n_seeds
  result <- cpp_order_hill_coverage_parallel(x, ordering$orders - 1L, q, n_cores)

  if (progress) cli_success("Done")

  hills <- result$hills
  coverage_mat <- result$coverage

  # Standardize at target coverage if requested
  standardized <- NULL
  if (!is.null(target_coverage)) {
    targets <- as.numeric(target_coverage)
    standardized <- lapply(seq_along(q), function(qi) {
      mat <- matrix(NA, nrow = n_seeds, ncol = length(targets))
      colnames(mat) <- paste0("C", targets)
      for (s in seq_len(n_seeds)) {
        mat[s, ] <- interpolate_at_coverage(hills[[qi]][s, ], coverage_mat[s, ], targets)
      }
      mat
    })
    names(standardized) <- names(hills)
  }

  structure(
    list(
      hills = hills,
      coverage = coverage_mat,
      pooled = colSums(x),
      standardized = standardized,
      q = q,
      target_coverage = target_coverage,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      distance = distance,
      method = method,
      focal_points = ordering$focal_points,
      call = match.call()
    ),
    class = "spacc_hill_coverage"
  )
}


# S3 methods ------------------------------------------------------------------

#' @export
print.spacc_hill_coverage <- function(x, ...) {
  .print_curves(x, "spacc Hill + Coverage",
                sprintf("Orders (q): %s", paste(x$q, collapse = ", ")),
                sprintf("Mean final coverage: %.3f", mean(x$coverage[, x$n_sites])),
                if (!is.null(x$target_coverage)) {
                  sprintf("Target coverage: %s", paste(x$target_coverage, collapse = ", "))
                })
}


#' @export
summary.spacc_hill_coverage <- function(object, ci_level = 0.95, ...) {
  cov_stats <- .curve_stats(object$coverage, ci_level)
  do.call(rbind, lapply(seq_along(object$q), function(qi) {
    hill_stats <- .curve_stats(object$hills[[qi]], ci_level)
    data.frame(
      q = object$q[qi],
      sites = seq_len(object$n_sites),
      mean_hill = hill_stats$mean,
      hill_lower = hill_stats$lower,
      hill_upper = hill_stats$upper,
      mean_coverage = cov_stats$mean,
      coverage_lower = cov_stats$lower,
      coverage_upper = cov_stats$upper
    )
  }))
}


#' @export
as.data.frame.spacc_hill_coverage <- function(x, row.names = NULL, optional = FALSE, ...) {
  summary(x)
}


#' @export
plot.spacc_hill_coverage <- function(x, xaxis = c("coverage", "sites"),
                                      ci = TRUE, ci_alpha = 0.2, ...) {
  xaxis <- match.arg(xaxis)
  summ <- summary(x)
  summ$q_label <- .q_label(summ$q)

  by_coverage <- xaxis == "coverage"
  .plot_curves_by(summ, "q_label", "Order",
                  x = if (by_coverage) "mean_coverage" else "sites",
                  y = "mean_hill", lower = "hill_lower", upper = "hill_upper",
                  xlab = if (by_coverage) "Sample coverage" else "Sites accumulated",
                  ylab = "Hill number",
                  title = if (by_coverage) "Hill Numbers vs. Sample Coverage"
                          else "Hill Numbers with Coverage Tracking",
                  ci = ci, ci_alpha = ci_alpha)
}
