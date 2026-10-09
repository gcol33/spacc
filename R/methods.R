# Helper to detect grouped spacc objects
is_grouped <- function(x) {
  is.list(x$curves) && !is.matrix(x$curves)
}


#' @export
print.spacc <- function(x, ...) {
  if (is_grouped(x)) {
    cat(sprintf("spacc: %d groups (%s)\n", length(x$group_names), x$method))
    for (i in seq_along(x$group_names)) {
      # Groups may have different site counts: read them off the curves
      cat(sprintf("  %s: %d sites, %d species\n",
                  x$group_names[i], ncol(x$curves[[i]]), x$n_species[[i]]))
    }
  } else {
    cat(sprintf(
      "spacc: %d sites, %d species, %d seeds (%s)\n",
      x$n_sites, x$n_species, x$n_seeds, x$method
    ))
  }
  invisible(x)
}


# Apply `fun(ungrouped_spacc, ...)` to every group of a grouped spacc object.
.lapply_groups <- function(object, fun, ...) {
  out <- lapply(seq_along(object$group_names), function(i) {
    tmp <- object
    tmp$curves <- object$curves[[i]]
    tmp$n_species <- object$n_species[[i]]
    tmp$n_sites <- ncol(tmp$curves)
    class(tmp) <- "spacc"
    fun(tmp, ...)
  })
  names(out) <- object$group_names
  out
}


#' @export
summary.spacc <- function(object, saturation_threshold = 0.9, ci_level = 0.95, ...) {

  if (is_grouped(object)) {
    return(.lapply_groups(object, summary, saturation_threshold = saturation_threshold,
                          ci_level = ci_level, ...))
  }

  n_sites <- object$n_sites
  stats <- .curve_stats(object$curves, ci_level, sd = TRUE)

  # Saturation point: first site where the mean reaches threshold * final mean
  threshold_value <- saturation_threshold * stats$mean[n_sites]
  saturation_point <- which(stats$mean >= threshold_value)[1]

  midpoint <- ceiling(n_sites / 2)
  cv_midpoint <- stats$sd[midpoint] / stats$mean[midpoint] * 100

  structure(
    list(
      sites = seq_len(n_sites),
      mean = stats$mean,
      lower = stats$lower,
      upper = stats$upper,
      sd = stats$sd,
      saturation_point = saturation_point,
      saturation_threshold = saturation_threshold,
      cv_midpoint = cv_midpoint,
      ci_level = ci_level,
      n_seeds = object$n_seeds,
      n_species = object$n_species,
      method = object$method
    ),
    class = "summary.spacc"
  )
}


#' @export
print.summary.spacc <- function(x, ...) {
  cat("Spatial Species Accumulation Curve\n")
  cat(strrep("-", 36), "\n")
  cat("Method:          ", x$method, "\n")
  cat("Seeds:           ", x$n_seeds, "\n")
  cat("Sites:           ", length(x$sites), "\n")
  cat("Total species:   ", x$n_species, "\n")
  cat(sprintf("Final species:   %.1f (%.0f%% CI: %.1f - %.1f)\n",
              x$mean[length(x$mean)],
              x$ci_level * 100,
              x$lower[length(x$lower)],
              x$upper[length(x$upper)]))
  cat(sprintf("Saturation (%d%%): %d sites\n",
              round(x$saturation_threshold * 100),
              x$saturation_point))
  cat(sprintf("CV at midpoint:  %.1f%%\n", x$cv_midpoint))
  invisible(x)
}


# Data frame of the curve statistics in a summary.spacc object.
.summary_frame <- function(s, group = NULL) {
  df <- data.frame(sites = s$sites, mean = s$mean, lower = s$lower,
                   upper = s$upper, sd = s$sd)
  if (!is.null(group)) df$group <- group
  df
}


#' @export
as.data.frame.spacc <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (is_grouped(x)) {
    summaries <- summary(x)
    return(do.call(rbind, lapply(names(summaries), function(g) {
      .summary_frame(summaries[[g]], group = g)
    })))
  }
  .summary_frame(summary(x))
}


#' Predict from Empirical Accumulation Curve
#'
#' Interpolate the mean empirical accumulation curve at arbitrary site counts
#' using linear interpolation. Unlike the predict method for `spacc_fit`
#' objects, this does not use a fitted model; it interpolates the observed
#' curve directly.
#'
#' @param object A `spacc` object.
#' @param n Numeric vector of site counts at which to interpolate. Defaults to
#'   25%, 50%, and 100% of total sites.
#' @param ci Logical. If `TRUE` (default), return a data frame with columns
#'   `n`, `mean`, `lower`, `upper`. If `FALSE`, return a named numeric vector.
#' @param ci_level Confidence level for the interval (default 0.95).
#' @param warn Logical. Warn when `n` values fall outside the observed range
#'   (default `TRUE`).
#' @param ... Ignored.
#'
#' @return A data frame (if `ci = TRUE`) or named numeric vector (if `ci = FALSE`).
#'   Out-of-range values return `NA`.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' sac <- spacc(species, coords, n_seeds = 10, progress = FALSE)
#' predict(sac, n = c(10, 25, 50))
#' predict(sac, n = c(10, 25), ci = FALSE)
#' }
#'
#' @export
predict.spacc <- function(object, n = NULL, ci = TRUE,
                          ci_level = 0.95, warn = TRUE, ...) {

  if (is_grouped(object)) {
    summaries <- summary(object, ci_level = ci_level)
    results <- lapply(names(summaries), function(g) {
      s <- summaries[[g]]
      res <- .predict_from_summary(s, n %||% .default_predict_sites(s), ci = ci, warn = warn)
      if (ci) {
        res$group <- g
      } else {
        names(res) <- paste0(g, ".", names(res))
      }
      res
    })
    return(if (ci) do.call(rbind, results) else unlist(results))
  }

  summ <- summary(object, ci_level = ci_level)
  .predict_from_summary(summ, n %||% .default_predict_sites(summ), ci = ci, warn = warn)
}


# Default interpolation points: a quarter, a half and all of the sites.
.default_predict_sites <- function(summ) {
  ns <- length(summ$sites)
  unique(c(floor(ns / 4), floor(ns / 2), ns))
}


# Internal interpolation helper
.predict_from_summary <- function(summ, n, ci, warn) {
  sites <- summ$sites

  out_of_range <- n < min(sites) | n > max(sites)
  if (warn && any(out_of_range)) {
    warning(sprintf(
      "n values outside observed range [%d, %d] will return NA: %s",
      min(sites), max(sites),
      paste(n[out_of_range], collapse = ", ")
    ))
  }

  interp_mean <- stats::approx(sites, summ$mean, xout = n, rule = 1)$y

  if (!ci) {
    out <- interp_mean
    names(out) <- n
    return(out)
  }

  data.frame(
    n = n,
    mean = interp_mean,
    lower = stats::approx(sites, summ$lower, xout = n, rule = 1)$y,
    upper = stats::approx(sites, summ$upper, xout = n, rule = 1)$y
  )
}


#' @export
as.data.frame.spacc_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  # x$data has columns x (sites) and y (observed mean)
  data.frame(
    sites = x$data$x,
    observed = x$data$y,
    predicted = unname(stats::predict(x$fit, newdata = data.frame(x = x$data$x))),
    asymptote = unname(x$asymptote),
    model = x$model
  )
}


#' @export
as.data.frame.spacc_hill <- function(x, row.names = NULL, optional = FALSE, ...) {
  .curve_frame(x$curves, "q", x$q, sd = TRUE)
}


#' @export
as.data.frame.spacc_phylo <- function(x, row.names = NULL, optional = FALSE, ...) {
  .curve_frame(x$curves, "metric", x$metric, sd = TRUE)
}


#' @export
as.data.frame.spacc_func <- function(x, row.names = NULL, optional = FALSE, ...) {
  .curve_frame(x$curves, "metric", x$metric, sd = TRUE)
}


#' @export
as.data.frame.spacc_dar <- function(x, row.names = NULL, optional = FALSE, ...) {
  summary(x, sd = TRUE)
}


# Column means and sds of the seed matrices named in `fields`.
.field_means <- function(x, fields, sd_fields = fields) {
  out <- data.frame(sites = seq_len(ncol(x[[fields[1]]])))
  for (f in fields) out[[f]] <- colMeans(x[[f]])
  for (f in sd_fields) out[[paste0(f, "_sd")]] <- apply(x[[f]], 2, stats::sd)
  out
}


#' @export
as.data.frame.spacc_beta <- function(x, row.names = NULL, optional = FALSE, ...) {
  .field_means(x, c("beta_total", "beta_turnover", "beta_nestedness"))
}


#' @export
as.data.frame.spacc_coverage <- function(x, row.names = NULL, optional = FALSE, ...) {
  .field_means(x, c("richness", "individuals", "coverage"),
               sd_fields = c("richness", "coverage"))
}


#' @export
as.data.frame.spacc_wavefront <- function(x, row.names = NULL, optional = FALSE, ...) {
  cbind(
    data.frame(step = seq_along(x$radius), radius = x$radius,
               sites = colMeans(x$sites_included)),
    .curve_stats(x$curves, sd = TRUE)
  )
}


#' @export
as.data.frame.spacc_decay <- function(x, row.names = NULL, optional = FALSE, ...) {
  cbind(data.frame(distance = x$breaks), .curve_stats(x$curves, sd = TRUE))
}


#' @export
as.data.frame.spacc_comp <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    comparison = paste(x$x_name, "vs", x$y_name),
    auc_x = x$auc_x,
    auc_y = x$auc_y,
    auc_diff = x$auc_diff,
    saturation_x = x$saturation_x,
    saturation_y = x$saturation_y,
    saturation_diff = x$saturation_diff,
    p_value = x$p_value,
    method = x$method
  )
}


#' @export
as.data.frame.spacc_metrics <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- x$metrics
  if (!is.null(x$coords)) {
    df$x <- x$coords$x
    df$y <- x$coords$y
  }
  df
}


#' @export
as.data.frame.spacc_endemism <- function(x, row.names = NULL, optional = FALSE, ...) {
  stats <- .curve_stats(x$endemism, sd = TRUE)
  data.frame(
    sites = seq_len(x$n_sites),
    mean_richness = colMeans(x$richness),
    mean_endemism = stats$mean,
    endemism_lower = stats$lower,
    endemism_upper = stats$upper,
    endemism_sd = stats$sd
  )
}


#' @export
as.data.frame.spacc_partition <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    q = x$q,
    alpha = x$alpha,
    beta = x$beta,
    gamma = x$gamma
  )
}


#' @export
as.data.frame.spacc_alpha <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$values
}


#' @export
`[.spacc` <- function(x, i, ...) {
  if (is_grouped(x)) {
    x$curves <- lapply(x$curves, function(m) m[i, , drop = FALSE])
    x$n_seeds <- nrow(x$curves[[1]])
  } else {
    x$curves <- x$curves[i, , drop = FALSE]
    if (!is.null(x$orders)) x$orders <- x$orders[i, , drop = FALSE]
    x$n_seeds <- nrow(x$curves)
  }
  x
}


#' Combine spacc Objects
#'
#' @param ... Named `spacc` objects to combine into a grouped `spacc`.
#' @return A grouped `spacc` object with per-group curves.
#' @export
c.spacc <- function(...) {
  objects <- list(...)
  nms <- names(objects)
  if (is.null(nms)) nms <- paste0("group_", seq_along(objects))
  .grouped_spacc(objects, nms, call = match.call())
}


#' @export
summary.spacc_decay <- function(object, ...) {
  df <- as.data.frame(object)
  structure(
    list(
      n_bins = length(object$breaks) - 1,
      n_seeds = object$n_seeds,
      distance_range = range(object$breaks),
      mean_similarity = df$mean,
      distance = df$distance
    ),
    class = "summary.spacc_decay"
  )
}


#' @export
print.summary.spacc_decay <- function(x, ...) {
  cat("Distance-Decay Summary\n")
  cat(strrep("-", 24), "\n")
  cat("Distance bins: ", x$n_bins, "\n")
  cat("Seeds:         ", x$n_seeds, "\n")
  cat(sprintf("Distance range: %.2f - %.2f\n", x$distance_range[1], x$distance_range[2]))
  cat(sprintf("Similarity range: %.3f - %.3f\n", min(x$mean_similarity), max(x$mean_similarity)))
  invisible(x)
}


#' @export
summary.spacc_rare <- function(object, ...) {
  structure(
    list(
      n_individuals = object$n,
      expected_richness = object$expected,
      se = object$sd,
      n_boot = object$n_boot
    ),
    class = "summary.spacc_rare"
  )
}


#' @export
print.summary.spacc_rare <- function(x, ...) {
  cat("Rarefaction Summary\n")
  cat(strrep("-", 20), "\n")
  n <- length(x$n_individuals)
  cat("Sample sizes: ", n, "\n")
  if (n > 0) {
    cat(sprintf("Range: %.0f - %.0f individuals\n", min(x$n_individuals), max(x$n_individuals)))
    cat(sprintf("Expected richness: %.1f - %.1f species\n",
                min(x$expected_richness), max(x$expected_richness)))
  }
  invisible(x)
}


#' Combine spacc_hill Objects
#'
#' Combine multiple `spacc_hill` objects by stacking their curve matrices.
#' All objects must have the same number of sites.
#'
#' @param ... `spacc_hill` objects to combine.
#' @return A combined `spacc_hill` object with more seeds.
#' @export
c.spacc_hill <- function(...) {
  .stack_seeds(list(...), "curves")
}


#' Combine spacc_beta Objects
#'
#' Combine multiple `spacc_beta` objects by stacking their curve matrices.
#' All objects must have the same number of sites and index.
#'
#' @param ... `spacc_beta` objects to combine.
#' @return A combined `spacc_beta` object with more seeds.
#' @export
c.spacc_beta <- function(...) {
  .stack_seeds(list(...), c("beta_total", "beta_turnover", "beta_nestedness", "distance"))
}


#' Combine spacc_coverage Objects
#'
#' Combine multiple `spacc_coverage` objects by stacking their curve matrices.
#' All objects must have the same number of sites.
#'
#' @param ... `spacc_coverage` objects to combine.
#' @return A combined `spacc_coverage` object with more seeds.
#' @export
c.spacc_coverage <- function(...) {
  .stack_seeds(list(...), c("richness", "individuals", "coverage"))
}
