# Shared helpers for the accumulation-curve classes.
#
# A curve object stores one or more seeds x steps matrices. These helpers turn
# such matrices into per-step summaries, long data frames, prints, plots and
# sf objects, so each class only states what is specific to it (which
# matrices, which key column, which labels).


# Resolve the `coords` argument of a front door: a data frame with x and y, an
# sf point layer, or a spacc_dist object carrying its own coordinates and
# distance matrix. Returns list(coords, dist_mat, crs).
.resolve_coords <- function(coords) {
  if (inherits(coords, "spacc_dist")) {
    return(list(coords = attr(coords, "coords"), dist_mat = as.matrix(coords), crs = NULL))
  }
  if (inherits(coords, "sf")) {
    check_suggests("sf")
    mat <- sf::st_coordinates(coords)
    return(list(coords = data.frame(x = mat[, 1], y = mat[, 2]), dist_mat = NULL,
                crs = sf::st_crs(coords)))
  }
  stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
  list(coords = coords, dist_mat = NULL, crs = NULL)
}


# Per-step mean, quantile band and (optionally) sd of a seeds x steps matrix.
.curve_stats <- function(mat, ci_level = 0.95, sd = FALSE) {
  alpha <- (1 - ci_level) / 2
  out <- data.frame(
    mean = colMeans(mat),
    lower = apply(mat, 2, stats::quantile, probs = alpha),
    upper = apply(mat, 2, stats::quantile, probs = 1 - alpha)
  )
  if (sd) out$sd <- apply(mat, 2, stats::sd)
  rownames(out) <- NULL
  out
}


# Long data frame over a list of curve matrices, one block per element, keyed
# by `key` = `values[i]`: columns key, sites, mean, lower, upper[, sd].
.curve_frame <- function(curves, key, values, ci_level = 0.95, sd = FALSE) {
  blocks <- lapply(seq_along(curves), function(i) {
    stats <- .curve_stats(curves[[i]], ci_level, sd)
    block <- data.frame(values[i], seq_len(nrow(stats)))
    names(block) <- c(key, "sites")
    cbind(block, stats)
  })
  do.call(rbind, blocks)
}


# Header line shared by the curve classes, followed by any extra lines.
.print_curves <- function(x, label, ...) {
  cat(sprintf("%s: %d sites, %d species, %d seeds\n",
              label, x$n_sites, x$n_species, x$n_seeds))
  for (line in c(...)) cat(line, "\n", sep = "")
  invisible(x)
}


# Factor of "q = <value>" labels in increasing q.
.q_label <- function(q) {
  factor(paste0("q = ", q), levels = paste0("q = ", sort(unique(q))))
}


# One line per level of `group`, with a quantile ribbon per level.
.plot_curves_by <- function(summ, group, group_label, xlab, ylab, title,
                            ci = TRUE, ci_alpha = 0.2, x = "sites", y = "mean",
                            lower = "lower", upper = "upper") {
  check_suggests("ggplot2")
  p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data[[x]], y = .data[[y]],
                                          color = .data[[group]]))
  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[lower]], ymax = .data[[upper]], fill = .data[[group]]),
      alpha = ci_alpha, color = NA
    )
  }
  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(x = xlab, y = ylab, color = group_label, fill = group_label,
                  title = title) +
    spacc_theme() +
    ggplot2::theme(legend.position = "bottom")
}


# Single mean curve with a quantile ribbon.
.plot_curve <- function(summ, xlab, ylab, title, subtitle = NULL,
                        ci = TRUE, ci_alpha = 0.3, x = "sites", y = "mean",
                        lower = "lower", upper = "upper") {
  check_suggests("ggplot2")
  p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[lower]], ymax = .data[[upper]]),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }
  p +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    spacc_theme()
}


# Constructor that produces each class with per-site map values (map = TRUE).
.map_constructor <- c(
  spacc_hill = "spaccHill",
  spacc_phylo = "spaccPhylo",
  spacc_func = "spaccFunc",
  spacc_beta = "spaccBeta",
  spacc_coverage = "spaccCoverage",
  spacc_endemism = "spaccEndemism"
)

.site_values_or_stop <- function(x) {
  if (is.null(x$site_values)) {
    constructor <- .map_constructor[class(x)[1]]
    stop(sprintf("No map data. Rerun %s() with map = TRUE.", constructor), call. = FALSE)
  }
  x$site_values
}

# Spatial map of one per-site column of `x$site_values`.
.plot_site_map <- function(x, column, title, subtitle, point_size, palette) {
  site_values <- .site_values_or_stop(x)
  if (!column %in% names(site_values)) {
    stop(sprintf("'%s' not found in site_values. Available: %s", column,
                 paste(setdiff(names(site_values), c("site_id", "x", "y")), collapse = ", ")),
         call. = FALSE)
  }
  plot_spatial_map(site_values, column, title = title, subtitle = subtitle,
                   point_size = point_size, palette = palette)
}

# as_sf() for every class whose per-site values live in `x$site_values`.
.as_sf_site_values <- function(x, crs = NULL) {
  as_sf_from_df(.site_values_or_stop(x), crs = crs)
}


# Stack the seed matrices named in `fields` across objects (lists of matrices
# are stacked element-wise) and add up the seeds. Used by the c() methods.
.stack_seeds <- function(objects, fields) {
  base <- objects[[1]]
  for (obj in objects[-1]) {
    if (obj$n_sites != base$n_sites) {
      stop("All objects must have the same number of sites")
    }
  }
  for (field in fields) {
    if (is.list(base[[field]])) {
      stacked <- lapply(seq_along(base[[field]]), function(i) {
        do.call(rbind, lapply(objects, function(obj) obj[[field]][[i]]))
      })
      names(stacked) <- names(base[[field]])
      base[[field]] <- stacked
    } else {
      base[[field]] <- do.call(rbind, lapply(objects, `[[`, field))
    }
  }
  base$n_seeds <- sum(vapply(objects, function(obj) obj$n_seeds, numeric(1)))
  base
}
