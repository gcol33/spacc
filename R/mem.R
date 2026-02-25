#' Spatial Eigenvectors (PCNM/dbMEM)
#'
#' Compute spatial eigenvectors from site coordinates using Principal
#' Coordinates of Neighbour Matrices (PCNM) or distance-based Moran's
#' Eigenvector Maps (dbMEM).
#'
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param threshold Numeric. Distance threshold for truncation. Default `NULL`
#'   uses the largest edge in the minimum spanning tree.
#' @param method Character. `"pcnm"` (default) retains all positive eigenvectors;
#'   `"dbmem"` retains only those with positive Moran's I.
#' @param distance Character. `"euclidean"` (default) or `"haversine"`.
#'
#' @return An object of class `spacc_mem` containing:
#'   \item{vectors}{Matrix of eigenvectors (sites x n_vectors)}
#'   \item{eigenvalues}{Positive eigenvalues}
#'   \item{moran_i}{Moran's I for each eigenvector}
#'   \item{threshold}{Truncation distance used}
#'   \item{coords}{Original coordinates}
#'   \item{n_sites}{Number of sites}
#'   \item{method}{Method used}
#'
#' @details
#' PCNM (Borcard & Legendre 2002) decomposes spatial structure into
#' orthogonal eigenvectors representing patterns at different spatial scales.
#' Large eigenvalues correspond to broad-scale patterns (positive spatial
#' autocorrelation), while small eigenvalues represent fine-scale patterns.
#'
#' The algorithm:
#' 1. Compute pairwise distances
#' 2. Truncate distances beyond threshold (set to 4 * threshold)
#' 3. Double-centre the squared truncated distance matrix
#' 4. Extract eigenvectors with positive eigenvalues
#' 5. For `"dbmem"`: additionally filter to Moran's I > 0
#'
#' @references
#' Borcard, D. & Legendre, P. (2002). All-scale spatial analysis of ecological
#' data by means of principal coordinates of neighbour matrices. Ecological
#' Modelling, 153, 51-68.
#'
#' Dray, S., Legendre, P. & Peres-Neto, P.R. (2006). Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbour
#' matrices (PCNM). Ecological Modelling, 196, 483-493.
#'
#' @seealso [spatialPartition()] for variance partitioning with MEMs
#'
#' @examples
#' coords <- data.frame(x = runif(30), y = runif(30))
#' mem <- spatialEigenvectors(coords)
#' print(mem)
#'
#' @export
spatialEigenvectors <- function(coords,
                                 threshold = NULL,
                                 method = c("pcnm", "dbmem"),
                                 distance = c("euclidean", "haversine")) {
  method <- match.arg(method)
  distance <- match.arg(distance)

  # Handle coords input

  if (inherits(coords, "spacc_dist")) {
    D <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    D <- as.matrix(cpp_distance_matrix(coords$x, coords$y, distance))
  }

  n_sites <- nrow(D)

  # Auto-threshold: largest edge in MST (Prim's algorithm)
  if (is.null(threshold)) {
    threshold <- mst_max_edge(D, n_sites)
  }

  # Truncate: D[D > threshold] <- 4 * threshold
  D_trunc <- D
  D_trunc[D_trunc > threshold] <- 4 * threshold

  # Double-centre the squared truncated distance matrix
  D2 <- D_trunc^2
  row_means <- rowMeans(D2)
  col_means <- colMeans(D2)
  grand_mean <- mean(D2)
  delta <- -0.5 * (sweep(sweep(D2, 1, row_means), 2, col_means) + grand_mean)

  # Eigen decomposition
  eig <- eigen(delta, symmetric = TRUE)

  # Retain positive eigenvalues only
  tol <- sqrt(.Machine$double.eps)
  pos <- eig$values > tol
  eigenvalues <- eig$values[pos]
  vectors <- eig$vectors[, pos, drop = FALSE]

  # Compute Moran's I for each eigenvector
  # Binary connectivity matrix: W[i,j] = 1 if D[i,j] <= threshold
  W <- (D <= threshold & D > 0) * 1
  W_sum <- sum(W)

  moran_i <- vapply(seq_len(ncol(vectors)), function(k) {
    v <- vectors[, k]
    v_centered <- v - mean(v)
    numerator <- n_sites / W_sum * sum(W * outer(v_centered, v_centered))
    denominator <- sum(v_centered^2)
    if (denominator == 0) return(0)
    numerator / denominator
  }, numeric(1))

  # For dbmem: filter to positive Moran's I only
  if (method == "dbmem") {
    keep <- moran_i > 0
    vectors <- vectors[, keep, drop = FALSE]
    eigenvalues <- eigenvalues[keep]
    moran_i <- moran_i[keep]
  }

  # Name the vectors
  colnames(vectors) <- paste0("MEM", seq_len(ncol(vectors)))

  structure(
    list(
      vectors = vectors,
      eigenvalues = eigenvalues,
      moran_i = moran_i,
      threshold = threshold,
      coords = coord_data,
      n_sites = n_sites,
      method = method,
      distance = distance
    ),
    class = "spacc_mem"
  )
}


# Prim's MST: return maximum edge weight
mst_max_edge <- function(D, n) {
  in_tree <- logical(n)
  in_tree[1] <- TRUE
  min_edge <- D[1, ]
  max_mst_edge <- 0

  for (step in seq_len(n - 1)) {
    # Find minimum edge to non-tree vertex
    best_j <- -1
    best_d <- Inf
    for (j in seq_len(n)) {
      if (!in_tree[j] && min_edge[j] < best_d) {
        best_d <- min_edge[j]
        best_j <- j
      }
    }

    in_tree[best_j] <- TRUE
    if (best_d > max_mst_edge) max_mst_edge <- best_d

    # Update minimum edges
    for (j in seq_len(n)) {
      if (!in_tree[j] && D[best_j, j] < min_edge[j]) {
        min_edge[j] <- D[best_j, j]
      }
    }
  }

  max_mst_edge
}


#' Spatial Variance Partitioning with MEMs
#'
#' Partition diversity variation into spatial and non-spatial components
#' using forward selection of Moran's Eigenvector Maps.
#'
#' @param x A spacc result object (e.g., `spacc_metrics`, `spacc_alpha`,
#'   `spacc_hill`) or a numeric vector of per-site diversity values.
#' @param mem A `spacc_mem` object from [spatialEigenvectors()].
#' @param metric Character. Which metric to extract from `x` (passed to
#'   internal extraction). Default `NULL` uses the first available.
#' @param forward Logical. Use forward selection? Default `TRUE`.
#'   If `FALSE`, all MEMs are included.
#' @param alpha Numeric. Significance threshold for forward selection.
#'   Default 0.05.
#'
#' @return An object of class `spacc_mem_partition` containing:
#'   \item{r_squared_spatial}{R-squared of the spatial model}
#'   \item{r_squared_total}{Total R-squared (same as spatial here)}
#'   \item{selected_mems}{Names of selected MEM vectors}
#'   \item{n_selected}{Number of selected MEMs}
#'   \item{anova_table}{ANOVA table from the final model}
#'   \item{coefficients}{Model coefficients}
#'
#' @details
#' Forward selection of MEMs proceeds by adding the MEM that most improves
#' the model AIC at each step, stopping when no MEM improves AIC by more
#' than 2 units or when p > alpha.
#'
#' @seealso [spatialEigenvectors()] for computing MEMs
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(30), y = runif(30))
#' species <- matrix(rpois(30 * 15, 2), nrow = 30)
#'
#' mem <- spatialEigenvectors(coords)
#' alpha <- alphaDiversity(species, q = 0)
#' part <- spatialPartition(alpha$q0, mem)
#' print(part)
#' }
#'
#' @export
spatialPartition <- function(x, mem, metric = NULL,
                              forward = TRUE, alpha = 0.05) {
  stopifnot(inherits(mem, "spacc_mem"))

  # Extract diversity values
  if (is.numeric(x) && is.null(dim(x))) {
    y <- x
  } else {
    y <- extract_partition_metric(x, metric)
  }

  stopifnot("x and mem must have same number of sites" = length(y) == mem$n_sites)

  vectors <- mem$vectors

  if (!forward || ncol(vectors) <= 1) {
    # Use all MEMs
    selected <- seq_len(ncol(vectors))
  } else {
    # Forward selection by AIC improvement
    selected <- integer(0)
    remaining <- seq_len(ncol(vectors))
    null_aic <- AIC(lm(y ~ 1))
    current_aic <- null_aic

    repeat {
      best_aic <- current_aic
      best_var <- NA

      for (v in remaining) {
        test_vars <- c(selected, v)
        test_df <- as.data.frame(vectors[, test_vars, drop = FALSE])
        test_fit <- lm(y ~ ., data = test_df)
        test_aic <- AIC(test_fit)

        if (test_aic < best_aic) {
          best_aic <- test_aic
          best_var <- v
        }
      }

      # Stop if no improvement (delta AIC > -2)
      if (is.na(best_var) || (current_aic - best_aic) < 2) break

      # Check significance
      test_vars <- c(selected, best_var)
      test_df <- as.data.frame(vectors[, test_vars, drop = FALSE])
      test_fit <- lm(y ~ ., data = test_df)
      pvals <- summary(test_fit)$coefficients[, 4]
      var_name <- colnames(vectors)[best_var]
      if (var_name %in% names(pvals) && pvals[var_name] > alpha) break

      selected <- c(selected, best_var)
      remaining <- setdiff(remaining, best_var)
      current_aic <- best_aic

      if (length(remaining) == 0) break
    }
  }

  if (length(selected) == 0) {
    return(structure(
      list(
        r_squared_spatial = 0,
        r_squared_total = 0,
        selected_mems = character(0),
        n_selected = 0,
        anova_table = NULL,
        coefficients = NULL
      ),
      class = "spacc_mem_partition"
    ))
  }

  # Final model
  final_df <- as.data.frame(vectors[, selected, drop = FALSE])
  final_fit <- lm(y ~ ., data = final_df)
  r2 <- summary(final_fit)$r.squared

  structure(
    list(
      r_squared_spatial = r2,
      r_squared_total = r2,
      selected_mems = colnames(vectors)[selected],
      n_selected = length(selected),
      anova_table = stats::anova(final_fit),
      coefficients = stats::coef(final_fit)
    ),
    class = "spacc_mem_partition"
  )
}


# Internal: extract per-site metric vector from various spacc objects
extract_partition_metric <- function(x, metric = NULL) {
  cls <- class(x)[1]

  switch(cls,
    "spacc_metrics" = {
      m_col <- if (!is.null(metric)) metric else x$metric_names[1]
      x$metrics[[m_col]]
    },
    "spacc_alpha" = {
      q_col <- if (!is.null(metric)) metric else paste0("q", x$q[1])
      x$values[[q_col]]
    },
    "spacc_hill" = {
      q_idx <- if (!is.null(metric)) match(metric, paste0("q", x$q)) else 1
      if (is.na(q_idx)) q_idx <- 1
      if (!is.null(x$site_values)) {
        x$site_values[[paste0("q", x$q[q_idx])]]
      } else {
        # Mean final value across seeds
        mat <- x$curves[[q_idx]]
        mat[, ncol(mat)]
      }
    },
    "spacc_partition" = {
      q_col <- if (!is.null(metric)) metric else paste0("q", x$q[1])
      x$alpha_per_site[[q_col]]
    },
    "data.frame" = {
      if (!is.null(metric) && metric %in% names(x)) {
        x[[metric]]
      } else {
        x[[1]]
      }
    },
    stop(sprintf("Cannot extract metric from class '%s'. Pass a numeric vector instead.", cls))
  )
}


# S3 methods for spacc_mem ---------------------------------------------------

#' @export
print.spacc_mem <- function(x, ...) {
  cat(sprintf("Spatial eigenvectors (%s): %d sites\n", toupper(x$method), x$n_sites))
  cat(sprintf("Threshold: %.4f\n", x$threshold))
  cat(sprintf("Eigenvectors retained: %d\n", ncol(x$vectors)))
  if (length(x$moran_i) > 0) {
    cat(sprintf("Moran's I range: [%.3f, %.3f]\n", min(x$moran_i), max(x$moran_i)))
  }
  invisible(x)
}


#' @export
summary.spacc_mem <- function(object, ...) {
  data.frame(
    vector = colnames(object$vectors),
    eigenvalue = object$eigenvalues,
    moran_i = object$moran_i,
    variance_explained = object$eigenvalues / sum(object$eigenvalues)
  )
}


#' @export
as.data.frame.spacc_mem <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- as.data.frame(x$vectors)
  if (!is.null(x$coords)) {
    df$x <- x$coords$x
    df$y <- x$coords$y
  }
  df
}


#' @export
plot.spacc_mem <- function(x, type = c("eigenvalues", "moran", "map"),
                            n_vectors = 4, point_size = 3,
                            palette = "viridis", ...) {
  type <- match.arg(type)
  check_suggests("ggplot2")

  if (type == "eigenvalues") {
    df <- data.frame(
      vector = seq_along(x$eigenvalues),
      eigenvalue = x$eigenvalues
    )
    ggplot2::ggplot(df, ggplot2::aes(x = vector, y = eigenvalue)) +
      ggplot2::geom_col(fill = "#4CAF50") +
      ggplot2::labs(x = "Eigenvector", y = "Eigenvalue",
                    title = sprintf("Spatial Eigenvalues (%s)", toupper(x$method))) +
      ggplot2::theme_minimal(base_size = 12)

  } else if (type == "moran") {
    df <- data.frame(
      vector = seq_along(x$moran_i),
      moran_i = x$moran_i
    )
    ggplot2::ggplot(df, ggplot2::aes(x = vector, y = moran_i)) +
      ggplot2::geom_col(fill = ifelse(df$moran_i > 0, "#4CAF50", "#F44336")) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "Eigenvector", y = "Moran's I",
                    title = "Moran's I per Eigenvector") +
      ggplot2::theme_minimal(base_size = 12)

  } else {
    # Map: show first n_vectors as spatial maps
    if (is.null(x$coords)) stop("No coordinates available.")

    n_show <- min(n_vectors, ncol(x$vectors))
    dfs <- lapply(seq_len(n_show), function(k) {
      data.frame(
        x = x$coords$x,
        y = x$coords$y,
        value = x$vectors[, k],
        vector = paste0("MEM", k)
      )
    })
    plot_df <- do.call(rbind, dfs)
    plot_df$vector <- factor(plot_df$vector, levels = paste0("MEM", seq_len(n_show)))

    ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = value)) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::scale_color_viridis_c(option = "C") +
      ggplot2::facet_wrap(~ vector) +
      ggplot2::labs(x = "x", y = "y", color = "Score",
                    title = "Spatial Eigenvector Maps") +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 12)
  }
}


# S3 methods for spacc_mem_partition ------------------------------------------

#' @export
print.spacc_mem_partition <- function(x, ...) {
  cat("Spatial Variance Partitioning\n")
  cat(strrep("-", 35), "\n")
  cat(sprintf("Spatial R-squared: %.3f\n", x$r_squared_spatial))
  cat(sprintf("Non-spatial:       %.3f\n", 1 - x$r_squared_spatial))
  cat(sprintf("MEMs selected:     %d\n", x$n_selected))
  if (x$n_selected > 0) {
    cat(sprintf("Selected: %s\n", paste(x$selected_mems, collapse = ", ")))
  }
  invisible(x)
}


#' @export
summary.spacc_mem_partition <- function(object, ...) {
  data.frame(
    component = c("Spatial", "Non-spatial"),
    r_squared = c(object$r_squared_spatial, 1 - object$r_squared_spatial),
    n_vectors = c(object$n_selected, NA)
  )
}


#' @export
plot.spacc_mem_partition <- function(x, ...) {
  check_suggests("ggplot2")

  df <- data.frame(
    component = c("Spatial", "Non-spatial"),
    r_squared = c(x$r_squared_spatial, 1 - x$r_squared_spatial)
  )
  df$component <- factor(df$component, levels = c("Spatial", "Non-spatial"))

  ggplot2::ggplot(df, ggplot2::aes(x = component, y = r_squared, fill = component)) +
    ggplot2::geom_col(width = 0.5) +
    ggplot2::scale_fill_manual(values = c("Spatial" = "#4CAF50", "Non-spatial" = "#78909C")) +
    ggplot2::labs(x = NULL, y = expression(R^2),
                  title = "Spatial Variance Partitioning") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")
}
