# Helper to detect grouped spacc objects
is_grouped <- function(x) {
  is.list(x$curves) && !is.matrix(x$curves)
}


#' @export
print.spacc <- function(x, ...) {
  if (is_grouped(x)) {
    cat(sprintf("spacc: %d groups (%s)\n", length(x$group_names), x$method))
    for (i in seq_along(x$group_names)) {
      cat(sprintf("  %s: %d sites, %d species\n",
                  x$group_names[i],
                  x$n_sites,
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
