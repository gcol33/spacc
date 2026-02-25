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
#'
#' @return An object of class `spacc_hill_coverage` containing:
#'   \item{hills}{Named list of n_seeds x n_sites matrices (one per q)}
#'   \item{coverage}{n_seeds x n_sites matrix of Chao-Jost coverage}
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
                               seed = NULL) {

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

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    if (progress) cli_info(sprintf("Computing distances (%d x %d)", nrow(x), nrow(x)))
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)

  storage.mode(x) <- "integer"

  if (progress) cli_info(sprintf("Computing Hill numbers + coverage (q = %s, %d seeds)",
                                  paste(q, collapse = ", "), n_seeds))

  result <- cpp_knn_hill_coverage_parallel(x, dist_mat, n_seeds, q, n_cores, progress)

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
        cov_vec <- coverage_mat[s, ]
        hill_vec <- hills[[qi]][s, ]
        mat[s, ] <- interpolate_at_coverage(hill_vec, cov_vec, targets)
      }
      mat
    })
    names(standardized) <- names(hills)
  }

  structure(
    list(
      hills = hills,
      coverage = coverage_mat,
      standardized = standardized,
      q = q,
      target_coverage = target_coverage,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_hill_coverage"
  )
}


# S3 methods ------------------------------------------------------------------

#' @export
print.spacc_hill_coverage <- function(x, ...) {
  mean_final_cov <- mean(x$coverage[, x$n_sites])
  cat(sprintf("spacc Hill + Coverage: %d sites, %d species, %d seeds\n",
              x$n_sites, x$n_species, x$n_seeds))
  cat(sprintf("Orders (q): %s\n", paste(x$q, collapse = ", ")))
  cat(sprintf("Mean final coverage: %.3f\n", mean_final_cov))
  if (!is.null(x$target_coverage)) {
    cat(sprintf("Target coverage: %s\n", paste(x$target_coverage, collapse = ", ")))
  }
  invisible(x)
}


#' @export
summary.spacc_hill_coverage <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2

  # Coverage summary
  cov_mean <- colMeans(object$coverage)
  cov_lower <- apply(object$coverage, 2, stats::quantile, probs = alpha)
  cov_upper <- apply(object$coverage, 2, stats::quantile, probs = 1 - alpha)

  results <- lapply(seq_along(object$q), function(qi) {
    mat <- object$hills[[qi]]
    data.frame(
      q = object$q[qi],
      sites = seq_len(object$n_sites),
      mean_hill = colMeans(mat),
      hill_lower = apply(mat, 2, stats::quantile, probs = alpha),
      hill_upper = apply(mat, 2, stats::quantile, probs = 1 - alpha),
      mean_coverage = cov_mean,
      coverage_lower = cov_lower,
      coverage_upper = cov_upper
    )
  })

  do.call(rbind, results)
}


#' @export
as.data.frame.spacc_hill_coverage <- function(x, row.names = NULL, optional = FALSE, ...) {
  summary(x)
}


#' @export
plot.spacc_hill_coverage <- function(x, xaxis = c("coverage", "sites"),
                                      ci = TRUE, ci_alpha = 0.2, ...) {
  xaxis <- match.arg(xaxis)
  check_suggests("ggplot2")

  summ <- summary(x)
  summ$q_label <- factor(paste0("q = ", summ$q),
                          levels = paste0("q = ", sort(unique(summ$q))))

  if (xaxis == "coverage") {
    p <- ggplot2::ggplot(summ, ggplot2::aes(x = mean_coverage, y = mean_hill,
                                             color = q_label))
    if (ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = hill_lower, ymax = hill_upper, fill = q_label),
        alpha = ci_alpha, color = NA
      )
    }
    p <- p +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(x = "Sample coverage", y = "Hill number",
                    color = "Order", fill = "Order",
                    title = "Hill Numbers vs. Sample Coverage")
  } else {
    p <- ggplot2::ggplot(summ, ggplot2::aes(x = sites, y = mean_hill,
                                             color = q_label))
    if (ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = hill_lower, ymax = hill_upper, fill = q_label),
        alpha = ci_alpha, color = NA
      )
    }
    p <- p +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(x = "Sites accumulated", y = "Hill number",
                    color = "Order", fill = "Order",
                    title = "Hill Numbers with Coverage Tracking")
  }

  p +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_hill_coverage)
autoplot.spacc_hill_coverage <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}
