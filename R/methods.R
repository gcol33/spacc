# Helper to detect grouped spacc objects
is_grouped <- function(x) {
  is.list(x$curves) && !is.matrix(x$curves)
}


#' @export
print.spacc <- function(x, ...) {
  if (is_grouped(x)) {
    cat(sprintf("spacc: %d groups (%s)\n", length(x$group_names), x$method))
    for (i in seq_along(x$group_names)) {
      # Get n_sites from curves matrix (groups may have different site counts)
      group_n_sites <- ncol(x$curves[[i]])
      cat(sprintf("  %s: %d sites, %d species\n",
                  x$group_names[i],
                  group_n_sites,
                  x$n_species[[i]]))
    }
  } else {
    cat(sprintf(
      "spacc: %d sites, %d species, %d seeds (%s)\n",
      x$n_sites, x$n_species, x$n_seeds, x$method
    ))
  }
  invisible(x)
}


#' @export
summary.spacc <- function(object, saturation_threshold = 0.9, ci_level = 0.95, ...) {

  if (is_grouped(object)) {
    # Summarize each group
    summaries <- lapply(seq_along(object$group_names), function(i) {
      # Create a temporary ungrouped spacc for summary
      tmp <- object
      tmp$curves <- object$curves[[i]]
      tmp$n_species <- object$n_species[[i]]
      # Get n_sites from the curves matrix (groups may have different site counts)
      tmp$n_sites <- ncol(tmp$curves)
      class(tmp) <- "spacc"
      summary(tmp, saturation_threshold = saturation_threshold, ci_level = ci_level, ...)
    })
    names(summaries) <- object$group_names
    return(summaries)
  }

  curves <- object$curves
  n_sites <- object$n_sites

  alpha <- (1 - ci_level) / 2

  # Compute statistics per site position
  mean_curve <- colMeans(curves)
  lower <- apply(curves, 2, stats::quantile, probs = alpha)
  upper <- apply(curves, 2, stats::quantile, probs = 1 - alpha)
  sd_curve <- apply(curves, 2, stats::sd)

  # Saturation point: first site where mean reaches threshold * max
  max_species <- mean_curve[n_sites]
  threshold_value <- saturation_threshold * max_species
  saturation_point <- which(mean_curve >= threshold_value)[1]

  # CV at midpoint
  midpoint <- ceiling(n_sites / 2)
  cv_midpoint <- sd_curve[midpoint] / mean_curve[midpoint] * 100

  structure(
    list(
      sites = seq_len(n_sites),
      mean = mean_curve,
      lower = lower,
      upper = upper,
      sd = sd_curve,
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


#' @export
as.data.frame.spacc <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (is_grouped(x)) {
    summaries <- summary(x)
    df_list <- lapply(names(summaries), function(g) {
      s <- summaries[[g]]
      data.frame(
        sites = s$sites,
        mean = s$mean,
        lower = s$lower,
        upper = s$upper,
        sd = s$sd,
        group = g
      )
    })
    return(do.call(rbind, df_list))
  }
  summ <- summary(x)
  data.frame(
    sites = summ$sites,
    mean = summ$mean,
    lower = summ$lower,
    upper = summ$upper,
    sd = summ$sd
  )
}


#' @export
as.data.frame.spacc_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  # x$data has columns x (sites) and y (observed mean)
  data.frame(
    sites = x$data$x,
    observed = x$data$y,
    predicted = stats::predict(x$fit, newdata = data.frame(x = x$data$x)),
    asymptote = x$asymptote,
    model = x$model
  )
}


#' @export
as.data.frame.spacc_hill <- function(x, row.names = NULL, optional = FALSE, ...) {
  # curves is a list of matrices, one per q value
  df_list <- lapply(seq_along(x$q), function(i) {
    curves <- x$curves[[i]]
    data.frame(
      sites = seq_len(ncol(curves)),
      q = x$q[i],
      mean = colMeans(curves),
      lower = apply(curves, 2, stats::quantile, probs = 0.025),
      upper = apply(curves, 2, stats::quantile, probs = 0.975),
      sd = apply(curves, 2, stats::sd)
    )
  })
  do.call(rbind, df_list)
}


#' @export
as.data.frame.spacc_beta <- function(x, row.names = NULL, optional = FALSE, ...) {
  # beta_total, beta_turnover, beta_nestedness are matrices (seeds x sites)
  data.frame(
    sites = seq_len(ncol(x$beta_total)),
    beta_total = colMeans(x$beta_total),
    beta_turnover = colMeans(x$beta_turnover),
    beta_nestedness = colMeans(x$beta_nestedness),
    beta_total_sd = apply(x$beta_total, 2, stats::sd),
    beta_turnover_sd = apply(x$beta_turnover, 2, stats::sd),
    beta_nestedness_sd = apply(x$beta_nestedness, 2, stats::sd)
  )
}


#' @export
as.data.frame.spacc_coverage <- function(x, row.names = NULL, optional = FALSE, ...) {
  # richness, individuals, coverage are matrices (seeds x sites)
  data.frame(
    sites = seq_len(ncol(x$richness)),
    richness = colMeans(x$richness),
    individuals = colMeans(x$individuals),
    coverage = colMeans(x$coverage),
    richness_sd = apply(x$richness, 2, stats::sd),
    coverage_sd = apply(x$coverage, 2, stats::sd)
  )
}


#' @export
as.data.frame.spacc_wavefront <- function(x, row.names = NULL, optional = FALSE, ...) {
  # curves is a matrix (seeds x steps), radius is a vector per step
  # sites_included is a matrix (seeds x steps) - take mean per step
  curves <- x$curves
  n_steps <- ncol(curves)

  # sites_included may be a matrix or need reshaping
  if (is.matrix(x$sites_included)) {
    sites_mean <- colMeans(x$sites_included)
  } else {
    # Reshape vector to matrix if needed
    sites_mean <- colMeans(matrix(x$sites_included, nrow = nrow(curves), ncol = n_steps))
  }

  data.frame(
    step = seq_len(n_steps),
    radius = x$radius,
    sites = sites_mean,
    mean = colMeans(curves),
    lower = apply(curves, 2, stats::quantile, probs = 0.025),
    upper = apply(curves, 2, stats::quantile, probs = 0.975),
    sd = apply(curves, 2, stats::sd)
  )
}


#' @export
as.data.frame.spacc_decay <- function(x, row.names = NULL, optional = FALSE, ...) {
  curves <- x$curves
  data.frame(
    distance = x$breaks,
    mean = colMeans(curves),
    lower = apply(curves, 2, stats::quantile, probs = 0.025),
    upper = apply(curves, 2, stats::quantile, probs = 0.975),
    sd = apply(curves, 2, stats::sd)
  )
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
as.data.frame.spacc_phylo <- function(x, row.names = NULL, optional = FALSE, ...) {
  curves <- x$curves
  data.frame(
    sites = seq_len(ncol(curves)),
    mean = colMeans(curves),
    lower = apply(curves, 2, stats::quantile, probs = 0.025),
    upper = apply(curves, 2, stats::quantile, probs = 0.975),
    sd = apply(curves, 2, stats::sd),
    metric = x$metric
  )
}


#' @export
as.data.frame.spacc_func <- function(x, row.names = NULL, optional = FALSE, ...) {
  curves <- x$curves
  data.frame(
    sites = seq_len(ncol(curves)),
    mean = colMeans(curves),
    lower = apply(curves, 2, stats::quantile, probs = 0.025),
    upper = apply(curves, 2, stats::quantile, probs = 0.975),
    sd = apply(curves, 2, stats::sd),
    metric = x$metric
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
as.data.frame.spacc_dar <- function(x, row.names = NULL, optional = FALSE, ...) {
  # DAR has Hill curves and area
  df_list <- lapply(seq_along(x$q), function(i) {
    data.frame(
      area = x$area,
      q = x$q[i],
      diversity = x$diversity[[i]]
    )
  })
  do.call(rbind, df_list)
}


#' @export
as.data.frame.spacc_endemism <- function(x, row.names = NULL, optional = FALSE, ...) {
  curves <- x$curves
  data.frame(
    sites = seq_len(ncol(curves)),
    mean = colMeans(curves),
    lower = apply(curves, 2, stats::quantile, probs = 0.025),
    upper = apply(curves, 2, stats::quantile, probs = 0.975),
    sd = apply(curves, 2, stats::sd),
    type = x$type
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
  df <- x$values
  df
}


#' @export
`[.spacc` <- function(x, i, ...) {
  if (is_grouped(x)) {
    x$curves <- lapply(x$curves, function(m) m[i, , drop = FALSE])
    x$n_seeds <- nrow(x$curves[[1]])
  } else {
    x$curves <- x$curves[i, , drop = FALSE]
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

  # Get names
  nms <- names(objects)
  if (is.null(nms)) {
    nms <- paste0("group_", seq_along(objects))
  }

  # Use first object as template
  base <- objects[[1]]

  structure(
    list(
      curves = stats::setNames(lapply(objects, `[[`, "curves"), nms),
      group_names = nms,
      coords = base$coords,
      n_seeds = base$n_seeds,
      n_sites = base$n_sites,
      n_species = stats::setNames(lapply(objects, `[[`, "n_species"), nms),
      method = base$method,
      distance = base$distance,
      backend = base$backend,
      sigma = base$sigma,
      cone_width = base$cone_width,
      time = base$time,
      w_space = base$w_space,
      w_time = base$w_time,
      support = base$support,
      call = match.call()
    ),
    class = "spacc"
  )
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
