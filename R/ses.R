#' Standardized Effect Size (SES) via Null Models
#'
#' Compute standardized effect sizes by comparing observed diversity metrics
#' against a null distribution from randomized communities. Supports multiple
#' null model algorithms and works with most spacc output classes.
#'
#' @param x A spacc output object. Supported classes: `spacc`, `spacc_hill`,
#'   `spacc_phylo`, `spacc_func`, `spacc_beta`, `spacc_metrics`,
#'   `spacc_alpha`, `spacc_partition`.
#' @param species A site-by-species matrix (required). The species matrix used
#'   to produce `x`. Needed because spacc objects do not store the raw species
#'   matrix.
#' @param coords Optional data.frame with columns `x` and `y`. Required if the
#'   original analysis used coordinates. If `x` contains stored coordinates,
#'   they will be used automatically.
#' @param metric Character or `NULL`. For multi-metric objects (e.g.,
#'   `spacc_hill` with multiple q), specify which metric to extract. If `NULL`,
#'   uses the first/default metric.
#' @param null_model Character. Null model algorithm:
#'   - `"frequency"`: Shuffle species columns independently (maintains column
#'     totals = species frequency)
#'   - `"richness"`: Shuffle species rows independently (maintains row totals
#'     = site richness)
#'   - `"both"`: Independent swap algorithm maintaining both row and column
#'     totals (Gotelli 2000)
#'   - `"curveball"`: Curveball algorithm for efficient swap (Strona et al. 2014)
#'   - `"torus"`: Toroidal shift preserving spatial autocorrelation. Shifts all
#'     coordinates by a random offset and reassigns species to shifted sites.
#'     Requires `coords`.
#'   - `"spatial_swap"`: Independent swap restricted to spatially proximate site
#'     pairs. Preserves both marginals while respecting spatial structure.
#'     Requires `coords`.
#' @param n_perm Integer. Number of permutations. Default 999.
#' @param parallel Logical. Use parallel processing for the underlying analysis?
#'   Default `TRUE`.
#' @param n_cores Integer. Number of cores. Default `NULL` (auto-detect).
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return An object of class `spacc_ses` containing:
#'   \item{observed}{Numeric vector of observed metric values}
#'   \item{null_mean}{Mean of null distribution}
#'   \item{null_sd}{Standard deviation of null distribution}
#'   \item{ses}{Standardized effect size: (observed - null_mean) / null_sd}
#'   \item{p_value}{Two-tailed p-value}
#'   \item{n_perm}{Number of permutations}
#'   \item{null_model}{Null model algorithm used}
#'   \item{metric}{Metric name}
#'   \item{input_class}{Class of input object}
#'
#' @details
#' SES is computed as:
#' \deqn{SES = \frac{observed - \bar{null}}{sd_{null}}}
#'
#' A two-tailed p-value is calculated as the proportion of null values at
#' least as extreme as the observed value:
#' \deqn{p = \frac{2 \cdot \min(r, n_{perm} + 1 - r)}{n_{perm} + 1}}
#' where \eqn{r} is the rank of the observed value among null values.
#'
#' **Null model algorithms:**
#' - `"frequency"`: Tests whether species composition matters given observed
#'   species frequencies
#' - `"richness"`: Tests whether species identity matters given observed site
#'   richness
#' - `"both"`: Maintains both marginal totals; tests non-random species
#'   co-occurrence patterns
#' - `"curveball"`: Efficient alternative to `"both"` with proven uniform
#'   sampling properties
#'
#' @references
#' Gotelli, N.J. (2000). Null model analysis of species co-occurrence patterns.
#' Ecology, 81, 2606-2621.
#'
#' Strona, G., Nappo, D., Boccacci, F., Fattorini, S. & San-Miguel-Ayanz, J.
#' (2014). A fast and unbiased procedure to randomize ecological binary
#' matrices with fixed row and column totals. Nature Communications, 5, 4114.
#'
#' @seealso [spaccHill()], [spaccBeta()], [spaccMetrics()]
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(20), y = runif(20))
#' species <- matrix(rbinom(20 * 15, 1, 0.3), nrow = 20)
#'
#' sac <- spacc(species, coords, n_seeds = 10)
#' result <- ses(sac, species, n_perm = 19)
#' print(result)
#' }
#'
#' @export
ses <- function(x, species,
                coords = NULL,
                metric = NULL,
                null_model = c("frequency", "richness", "both", "curveball",
                               "torus", "spatial_swap"),
                n_perm = 999L,
                parallel = TRUE,
                n_cores = NULL,
                progress = TRUE,
                seed = NULL) {

  null_model <- match.arg(null_model)
  if (!is.null(seed)) set.seed(seed)

  species <- as.matrix(species)
  input_class <- class(x)[1]

  # Spatial null models require coords
  spatial_nulls <- c("torus", "spatial_swap")
  if (null_model %in% spatial_nulls && is.null(coords) && is.null(x$coords)) {
    stop(sprintf("Null model '%s' requires coords.", null_model))
  }

  # Supported classes

  supported <- c("spacc", "spacc_hill", "spacc_phylo", "spacc_func",
                 "spacc_beta", "spacc_metrics", "spacc_alpha", "spacc_partition")
  if (!input_class %in% supported) {
    stop(sprintf("ses() does not support class '%s'. Supported: %s",
                 input_class, paste(supported, collapse = ", ")))
  }

  # Get coords from object if not provided
  if (is.null(coords) && !is.null(x$coords)) {
    coords <- x$coords
  }

  # Extract observed values
  observed <- extract_ses_metric(x, metric)
  metric_name <- attr(observed, "metric_name")

  if (progress) cli_info(sprintf("Running %d null model permutations (%s)", n_perm, null_model))

  # Run permutations
  null_values <- matrix(NA, nrow = n_perm, ncol = length(observed))

  for (p in seq_len(n_perm)) {
    # Randomize species matrix
    null_species <- randomize_matrix(species, null_model, coords = coords)

    # Re-run analysis with null species
    null_result <- tryCatch(
      rerun_analysis(x, null_species, coords, parallel, n_cores),
      error = function(e) NULL
    )

    if (!is.null(null_result)) {
      null_values[p, ] <- extract_ses_metric(null_result, metric)
    }
  }

  # Remove failed permutations
  valid <- !is.na(null_values[, 1])
  null_values <- null_values[valid, , drop = FALSE]
  n_valid <- nrow(null_values)

  if (n_valid < 5) {
    stop("Too few valid permutations. Check that the analysis can run on randomized data.")
  }

  # Compute SES
  null_mean <- colMeans(null_values)
  null_sd <- apply(null_values, 2, stats::sd)

  ses_vals <- ifelse(null_sd > 0, (observed - null_mean) / null_sd, 0)

  # Two-tailed p-value
  p_value <- vapply(seq_along(observed), function(i) {
    rank_obs <- sum(null_values[, i] <= observed[i]) + 1
    2 * min(rank_obs, n_valid + 2 - rank_obs) / (n_valid + 1)
  }, numeric(1))

  if (progress) cli_success(sprintf("Done (%d valid permutations)", n_valid))

  structure(
    list(
      observed = observed,
      null_mean = null_mean,
      null_sd = null_sd,
      ses = ses_vals,
      p_value = p_value,
      n_perm = n_perm,
      n_valid = n_valid,
      null_model = null_model,
      metric = metric_name,
      input_class = input_class,
      coords = coords,
      call = match.call()
    ),
    class = "spacc_ses"
  )
}


# Extract the relevant metric from a spacc result object
extract_ses_metric <- function(x, metric = NULL) {
  cls <- class(x)[1]

  result <- switch(cls,
    "spacc" = {
      # Mean curve across seeds
      vals <- colMeans(x$curves)
      attr(vals, "metric_name") <- "richness"
      vals
    },
    "spacc_hill" = {
      # Pick first q or specified metric
      q_idx <- if (!is.null(metric)) {
        match(metric, paste0("q", x$q))
      } else 1
      if (is.na(q_idx)) q_idx <- 1
      vals <- colMeans(x$curves[[q_idx]])
      attr(vals, "metric_name") <- paste0("hill_q", x$q[q_idx])
      vals
    },
    "spacc_phylo" = {
      m_idx <- if (!is.null(metric)) match(metric, x$metric) else 1
      if (is.na(m_idx)) m_idx <- 1
      vals <- colMeans(x$curves[[m_idx]])
      attr(vals, "metric_name") <- x$metric[m_idx]
      vals
    },
    "spacc_func" = {
      m_idx <- if (!is.null(metric)) match(metric, x$metric) else 1
      if (is.na(m_idx)) m_idx <- 1
      vals <- colMeans(x$curves[[m_idx]])
      attr(vals, "metric_name") <- x$metric[m_idx]
      vals
    },
    "spacc_beta" = {
      vals <- colMeans(x$beta_total)
      attr(vals, "metric_name") <- "beta_total"
      vals
    },
    "spacc_metrics" = {
      m_col <- if (!is.null(metric)) metric else x$metric_names[1]
      vals <- x$metrics[[m_col]]
      attr(vals, "metric_name") <- m_col
      vals
    },
    "spacc_alpha" = {
      q_col <- if (!is.null(metric)) metric else paste0("q", x$q[1])
      vals <- x$values[[q_col]]
      attr(vals, "metric_name") <- q_col
      vals
    },
    "spacc_partition" = {
      vals <- x$gamma
      attr(vals, "metric_name") <- "gamma"
      vals
    },
    stop(sprintf("Unsupported class: %s", cls))
  )

  result
}


# Re-run the same analysis with a new species matrix
rerun_analysis <- function(original, null_species, coords, parallel, n_cores) {
  cls <- class(original)[1]
  n_cores_val <- resolve_cores(n_cores, parallel)

  switch(cls,
    "spacc" = {
      spacc(null_species, coords,
            n_seeds = original$n_seeds,
            method = original$method,
            distance = original$distance,
            parallel = parallel, n_cores = n_cores,
            progress = FALSE)
    },
    "spacc_hill" = {
      spaccHill(null_species, coords,
                q = original$q,
                n_seeds = original$n_seeds,
                method = original$method,
                distance = original$distance,
                parallel = parallel, n_cores = n_cores,
                progress = FALSE)
    },
    "spacc_beta" = {
      spaccBeta(null_species, coords,
                n_seeds = original$n_seeds,
                method = original$method,
                index = original$index,
                parallel = parallel, n_cores = n_cores,
                progress = FALSE)
    },
    "spacc_alpha" = {
      alphaDiversity(null_species, q = original$q, coords = coords)
    },
    "spacc_partition" = {
      diversityPartition(null_species, q = original$q, coords = coords)
    },
    "spacc_metrics" = {
      spaccMetrics(null_species, coords,
                   metrics = original$metric_names,
                   method = original$method,
                   parallel = parallel, n_cores = n_cores,
                   progress = FALSE)
    },
    # For phylo/func we can't easily re-run without tree/traits
    # Return NULL to skip
    NULL
  )
}


# Null model algorithms --------------------------------------------------

randomize_matrix <- function(m, method, coords = NULL) {
  pa <- (m > 0) * 1L

  switch(method,
    "frequency" = {
      # Shuffle each column independently
      for (j in seq_len(ncol(pa))) {
        pa[, j] <- sample(pa[, j])
      }
      pa
    },
    "richness" = {
      # Shuffle each row independently
      for (i in seq_len(nrow(pa))) {
        pa[i, ] <- sample(pa[i, ])
      }
      pa
    },
    "both" = {
      # Independent swap: maintain both marginals
      independent_swap(pa, n_swaps = 5L * sum(pa))
    },
    "curveball" = {
      curveball(pa)
    },
    "torus" = {
      torus_shift(pa, coords)
    },
    "spatial_swap" = {
      spatial_swap(pa, coords)
    }
  )
}


# Independent swap algorithm (Gotelli 2000)
independent_swap <- function(m, n_swaps) {
  nr <- nrow(m)
  nc <- ncol(m)

  for (s in seq_len(n_swaps)) {
    # Pick 2 random rows and 2 random columns
    rows <- sample(nr, 2)
    cols <- sample(nc, 2)

    submat <- m[rows, cols]

    # Check for 2x2 checkerboard pattern to swap
    if ((submat[1, 1] == 1 && submat[2, 2] == 1 &&
         submat[1, 2] == 0 && submat[2, 1] == 0) ||
        (submat[1, 1] == 0 && submat[2, 2] == 0 &&
         submat[1, 2] == 1 && submat[2, 1] == 1)) {
      # Swap
      m[rows[1], cols[1]] <- 1L - m[rows[1], cols[1]]
      m[rows[2], cols[2]] <- 1L - m[rows[2], cols[2]]
      m[rows[1], cols[2]] <- 1L - m[rows[1], cols[2]]
      m[rows[2], cols[1]] <- 1L - m[rows[2], cols[1]]
    }
  }
  m
}


# Curveball algorithm (Strona et al. 2014)
curveball <- function(m) {
  nr <- nrow(m)

  # Get column indices of 1s for each row
  row_cols <- lapply(seq_len(nr), function(i) which(m[i, ] == 1L))

  n_trades <- 5L * nr

  for (t in seq_len(n_trades)) {
    # Pick two random rows
    ab <- sample(nr, 2)
    a <- ab[1]
    b <- ab[2]

    # Find exclusive columns
    only_a <- setdiff(row_cols[[a]], row_cols[[b]])
    only_b <- setdiff(row_cols[[b]], row_cols[[a]])

    n_a <- length(only_a)
    n_b <- length(only_b)
    n_trade <- min(n_a, n_b)
    if (n_trade == 0) next

    # Randomly select number and which columns to trade
    # safe_sample avoids sample(x, n) interpreting scalar x as 1:x
    n_swap <- sample.int(n_trade, 1)
    swap_a <- if (n_a == 1) only_a else only_a[sample.int(n_a, n_swap)]
    swap_b <- if (n_b == 1) only_b else only_b[sample.int(n_b, n_swap)]

    # Trade: row a gives swap_a, receives swap_b; row b gives swap_b, receives swap_a
    row_cols[[a]] <- sort(c(setdiff(row_cols[[a]], swap_a), swap_b))
    row_cols[[b]] <- sort(c(setdiff(row_cols[[b]], swap_b), swap_a))
  }

  # Reconstruct matrix
  result <- matrix(0L, nrow = nr, ncol = ncol(m))
  for (i in seq_len(nr)) {
    result[i, row_cols[[i]]] <- 1L
  }
  result
}


# Toroidal shift null model
# Shifts all coordinates by random offset, then reassigns species
# to shifted sites via nearest-neighbor mapping
torus_shift <- function(pa, coords) {
  n <- nrow(pa)
  x_range <- range(coords$x)
  y_range <- range(coords$y)
  x_width <- diff(x_range)
  y_width <- diff(y_range)

  # Random shift
  dx <- stats::runif(1, 0, x_width)
  dy <- stats::runif(1, 0, y_width)

  # Shift coordinates with wrapping
  new_x <- ((coords$x - x_range[1] + dx) %% x_width) + x_range[1]
  new_y <- ((coords$y - y_range[1] + dy) %% y_width) + y_range[1]

  # For each shifted site, find nearest original site
  # Simple O(n^2) matching - fine for typical use
  mapping <- integer(n)
  for (i in seq_len(n)) {
    dists <- (new_x[i] - coords$x)^2 + (new_y[i] - coords$y)^2
    mapping[i] <- which.min(dists)
  }

  # Reassign rows
  pa[mapping, , drop = FALSE]
}


# Spatially-constrained independent swap
# Like independent_swap but only swaps between spatially proximate sites
spatial_swap <- function(pa, coords, n_swaps = NULL) {
  nr <- nrow(pa)
  nc <- ncol(pa)

  # Compute pairwise distances
  dist_vec <- as.matrix(stats::dist(cbind(coords$x, coords$y)))

  # Distance threshold: median pairwise distance
  threshold <- stats::median(dist_vec[upper.tri(dist_vec)])

  # Find valid pairs within threshold
  valid_pairs <- which(dist_vec <= threshold & upper.tri(dist_vec), arr.ind = TRUE)
  n_pairs <- nrow(valid_pairs)

  if (n_pairs == 0) {
    # Fallback to standard independent swap
    return(independent_swap(pa, n_swaps = 5L * sum(pa)))
  }

  if (is.null(n_swaps)) n_swaps <- 5L * sum(pa)

  for (s in seq_len(n_swaps)) {
    # Pick a random proximate pair
    pair_idx <- sample.int(n_pairs, 1)
    rows <- valid_pairs[pair_idx, ]
    cols <- sample(nc, 2)

    submat <- pa[rows, cols]

    # Check for checkerboard pattern to swap
    if ((submat[1, 1] == 1 && submat[2, 2] == 1 &&
         submat[1, 2] == 0 && submat[2, 1] == 0) ||
        (submat[1, 1] == 0 && submat[2, 2] == 0 &&
         submat[1, 2] == 1 && submat[2, 1] == 1)) {
      pa[rows[1], cols[1]] <- 1L - pa[rows[1], cols[1]]
      pa[rows[2], cols[2]] <- 1L - pa[rows[2], cols[2]]
      pa[rows[1], cols[2]] <- 1L - pa[rows[1], cols[2]]
      pa[rows[2], cols[1]] <- 1L - pa[rows[2], cols[1]]
    }
  }
  pa
}


# S3 methods for spacc_ses -----------------------------------------------

#' @export
print.spacc_ses <- function(x, ...) {
  cat(sprintf("Standardized Effect Size (SES)\n"))
  cat(strrep("-", 36), "\n")
  cat(sprintf("Input class: %s\n", x$input_class))
  cat(sprintf("Metric: %s\n", x$metric))
  cat(sprintf("Null model: %s (%d permutations)\n", x$null_model, x$n_valid))
  cat(sprintf("SES range: [%.2f, %.2f]\n", min(x$ses), max(x$ses)))
  cat(sprintf("Mean SES: %.2f\n", mean(x$ses)))

  n_sig <- sum(x$p_value < 0.05)
  cat(sprintf("Significant (p < 0.05): %d / %d values (%.0f%%)\n",
              n_sig, length(x$p_value), 100 * n_sig / length(x$p_value)))
  invisible(x)
}


#' @export
summary.spacc_ses <- function(object, ...) {
  data.frame(
    position = seq_along(object$observed),
    observed = object$observed,
    null_mean = object$null_mean,
    null_sd = object$null_sd,
    ses = object$ses,
    p_value = object$p_value
  )
}


#' @export
as.data.frame.spacc_ses <- function(x, row.names = NULL, optional = FALSE, ...) {
  summary(x)
}


#' @export
as_sf.spacc_ses <- function(x, crs = NULL) {
  if (is.null(x$coords)) stop("No coordinates available. ses() input had no coords.")
  # Only meaningful for per-site input classes
  per_site_classes <- c("spacc_metrics", "spacc_alpha", "spacc_partition")
  if (!x$input_class %in% per_site_classes) {
    stop(sprintf("as_sf() only meaningful for per-site SES (input class: %s). Supported: %s",
                 x$input_class, paste(per_site_classes, collapse = ", ")))
  }
  df <- data.frame(
    x = x$coords$x,
    y = x$coords$y,
    ses = x$ses,
    p_value = x$p_value,
    observed = x$observed,
    null_mean = x$null_mean
  )
  as_sf_from_df(df, crs = crs)
}


#' @export
plot.spacc_ses <- function(x, type = c("curve", "histogram", "map"),
                            significance = 1.96,
                            metric_map = c("ses", "p_value"),
                            point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$coords)) stop("No coordinates available. ses() input had no coords.")
    metric_map <- match.arg(metric_map)
    vals <- if (metric_map == "ses") x$ses else x$p_value
    df <- data.frame(x = x$coords$x, y = x$coords$y)
    df[[metric_map]] <- vals
    return(plot_spatial_map(df, metric_map,
                            title = sprintf("SES %s: %s", metric_map, x$metric),
                            subtitle = sprintf("%s null model, %d permutations",
                                               x$null_model, x$n_valid),
                            point_size = point_size, palette = palette))
  }

  check_suggests("ggplot2")

  if (type == "curve") {
    df <- data.frame(
      position = seq_along(x$ses),
      ses = x$ses,
      p_value = x$p_value
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["position"]], y = .data[["ses"]])) +
      ggplot2::geom_hline(yintercept = c(-significance, 0, significance),
                           linetype = c("dashed", "solid", "dashed"),
                           color = c("#C62828", "grey50", "#C62828")) +
      ggplot2::geom_line(linewidth = 0.8, color = "#2E7D32") +
      ggplot2::geom_point(
        data = df[df$p_value < 0.05, ],
        ggplot2::aes(x = .data[["position"]], y = .data[["ses"]]),
        color = "#C62828", size = 1.5
      ) +
      ggplot2::labs(
        x = "Position",
        y = "Standardized Effect Size",
        title = sprintf("SES: %s (%s null model)", x$metric, x$null_model)
      ) +
      ggplot2::theme_minimal(base_size = 12)

    return(p)

  } else {
    # Histogram of SES values
    df <- data.frame(ses = x$ses)

    ggplot2::ggplot(df, ggplot2::aes(x = .data[["ses"]])) +
      ggplot2::geom_histogram(bins = 30, fill = "#4CAF50", color = "white", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = c(-significance, significance),
                           linetype = "dashed", color = "#C62828") +
      ggplot2::geom_vline(xintercept = 0, color = "grey50") +
      ggplot2::labs(
        x = "Standardized Effect Size",
        y = "Count",
        title = sprintf("SES Distribution: %s", x$metric)
      ) +
      ggplot2::theme_minimal(base_size = 12)
  }
}
