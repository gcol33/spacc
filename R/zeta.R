#' Zeta Diversity
#'
#' Compute zeta diversity — the mean number of species shared across k sites —
#' for increasing orders of k. The zeta decline curve reveals community assembly
#' processes: exponential decline suggests stochastic assembly, while power-law
#' decline indicates niche-based assembly.
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#'   Automatically binarized.
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param orders Integer vector. Orders of zeta diversity to compute (number of
#'   sites in each combination). Default `1:10`.
#' @param n_samples Integer. Number of random combinations to sample per order.
#'   Default 100.
#' @param method Character. Method for selecting k-site combinations:
#'   `"knn"` (spatially nearest sites) or `"random"` (random combinations).
#'   Default `"knn"`.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param seed Integer. Random seed for reproducibility. Default `NULL`.
#' @param progress Logical. Show progress? Default `TRUE`.
#'
#' @return An object of class `spacc_zeta` containing:
#'   \item{zeta}{Mean zeta values per order}
#'   \item{zeta_sd}{Standard deviations per order}
#'   \item{orders}{The k values}
#'   \item{n_samples}{Number of samples per order}
#'   \item{ratio}{Zeta ratio: zeta_k / zeta_(k-1)}
#'   \item{decline}{Data.frame with exponential and power-law fit statistics}
#'   \item{method}{Method used}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Total species count}
#'
#' @details
#' Zeta diversity of order k (\eqn{\zeta_k}) is the mean number of species
#' shared across k sites. Key properties:
#'
#' - \eqn{\zeta_1} = mean species richness per site
#' - \eqn{\zeta_2} = mean number of species shared by any two sites
#' - \eqn{\zeta_k} decreases monotonically with k
#'
#' The zeta decline ratio (\eqn{\zeta_k / \zeta_{k-1}}) is diagnostic:
#' - Constant ratio: exponential decline (stochastic assembly)
#' - Increasing ratio: power-law decline (deterministic/niche-based assembly)
#'
#' The `knn` method selects spatially nearest k sites from each focal site,
#' which is ecologically meaningful for testing spatial turnover. The `random`
#' method samples random k-site combinations, providing a null expectation.
#'
#' @references
#' Hui, C. & McGeoch, M.A. (2014). Zeta diversity as a concept and metric
#' that unifies incidence-based biodiversity patterns. The American Naturalist,
#' 184, 684-694.
#'
#' Latombe, G., McGeoch, M.A., Nipperess, D.A. & Hui, C. (2018). zetadiv:
#' an R package for computing compositional change across multiple sites,
#' assemblages or cases. bioRxiv, 324897.
#'
#' @seealso [spaccBeta()] for pairwise beta diversity, [distanceDecay()] for
#'   distance-decay relationships
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(30), y = runif(30))
#' species <- matrix(rbinom(30 * 20, 1, 0.3), nrow = 30)
#' zeta <- zetaDiversity(species, coords, orders = 1:5, n_samples = 50)
#' plot(zeta)
#' }
#'
#' @export
zetaDiversity <- function(x, coords,
                           orders = 1:10,
                           n_samples = 100L,
                           method = c("knn", "random"),
                           distance = c("euclidean", "haversine"),
                           seed = NULL,
                           progress = TRUE) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  x <- as.matrix(x)
  pa <- (x > 0) * 1L

  n_sites <- nrow(pa)
  n_species <- ncol(pa)

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    dist_mat <- NULL
  }

  stopifnot("x and coords must have same number of rows" = nrow(x) == nrow(coord_data))

  # Clip orders to available sites
  orders <- orders[orders >= 1 & orders <= n_sites]
  if (length(orders) == 0) stop("No valid orders. Must be between 1 and n_sites.")

  # Compute distance matrix if needed for knn
  if (method == "knn" && is.null(dist_mat)) {
    if (progress) cli_info("Computing distance matrix")
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  if (!is.null(seed)) set.seed(seed)

  if (progress) cli_info(sprintf("Computing zeta diversity (orders %d-%d, %s)",
                                  min(orders), max(orders), method))

  zeta_vals <- numeric(length(orders))
  zeta_sds <- numeric(length(orders))

  for (idx in seq_along(orders)) {
    k <- orders[idx]

    if (k == 1) {
      # Zeta_1 = mean species richness per site
      richness <- rowSums(pa)
      zeta_vals[idx] <- mean(richness)
      zeta_sds[idx] <- stats::sd(richness)
      next
    }

    # Generate k-site combinations and compute shared species
    shared_counts <- numeric(n_samples)

    for (s in seq_len(n_samples)) {
      if (method == "knn") {
        # Pick a random focal site, then its k-1 nearest neighbors
        focal <- sample.int(n_sites, 1)
        neighbors <- order(dist_mat[focal, ])[seq_len(k)]  # includes focal
        sites_idx <- neighbors
      } else {
        # Random k sites
        sites_idx <- sample.int(n_sites, k)
      }

      # Count species present in ALL k sites
      shared <- pa[sites_idx, , drop = FALSE]
      shared_counts[s] <- sum(colSums(shared) == k)
    }

    zeta_vals[idx] <- mean(shared_counts)
    zeta_sds[idx] <- stats::sd(shared_counts)
  }

  # Compute zeta ratio
  ratio <- rep(NA_real_, length(orders))
  for (i in 2:length(orders)) {
    if (zeta_vals[i - 1] > 0) {
      ratio[i] <- zeta_vals[i] / zeta_vals[i - 1]
    }
  }

  # Fit decline models (exponential and power-law)
  decline <- fit_zeta_decline(orders, zeta_vals)

  if (progress) cli_success("Done")

  structure(
    list(
      zeta = zeta_vals,
      zeta_sd = zeta_sds,
      orders = orders,
      n_samples = n_samples,
      ratio = ratio,
      decline = decline,
      method = method,
      distance = distance,
      coords = coord_data,
      n_sites = n_sites,
      n_species = n_species
    ),
    class = "spacc_zeta"
  )
}


#' Fit exponential and power-law decline models to zeta values
#' @noRd
fit_zeta_decline <- function(orders, zeta_vals) {
  # Need at least 3 non-zero points
  valid <- zeta_vals > 0
  if (sum(valid) < 3) {
    return(data.frame(model = c("exponential", "power_law"),
                      r_squared = c(NA, NA),
                      aic = c(NA, NA)))
  }

  k <- orders[valid]
  z <- zeta_vals[valid]

  # Exponential: log(zeta) ~ k => zeta = a * exp(b * k)
  exp_fit <- tryCatch({
    lm_fit <- stats::lm(log(z) ~ k)
    list(r_sq = summary(lm_fit)$r.squared,
         aic = stats::AIC(lm_fit))
  }, error = function(e) list(r_sq = NA, aic = NA))

  # Power-law: log(zeta) ~ log(k) => zeta = a * k^b
  pow_fit <- tryCatch({
    lm_fit <- stats::lm(log(z) ~ log(k))
    list(r_sq = summary(lm_fit)$r.squared,
         aic = stats::AIC(lm_fit))
  }, error = function(e) list(r_sq = NA, aic = NA))

  data.frame(
    model = c("exponential", "power_law"),
    r_squared = c(exp_fit$r_sq, pow_fit$r_sq),
    aic = c(exp_fit$aic, pow_fit$aic)
  )
}


# S3 methods ---------------------------------------------------------------

#' @export
print.spacc_zeta <- function(x, ...) {
  cat(sprintf("Zeta Diversity: %d sites, %d species\n", x$n_sites, x$n_species))
  cat(sprintf("Method: %s, %d samples per order\n", x$method, x$n_samples))
  cat(strrep("-", 40), "\n")

  df <- data.frame(
    order = x$orders,
    zeta = round(x$zeta, 2),
    sd = round(x$zeta_sd, 2),
    ratio = round(x$ratio, 3)
  )
  print(df, row.names = FALSE)

  if (!all(is.na(x$decline$aic))) {
    cat("\nDecline model fits:\n")
    print(x$decline, row.names = FALSE)
    best <- x$decline$model[which.min(x$decline$aic)]
    cat(sprintf("Best fit: %s\n", best))
  }

  invisible(x)
}


#' @export
summary.spacc_zeta <- function(object, ...) {
  data.frame(
    order = object$orders,
    zeta = object$zeta,
    sd = object$zeta_sd,
    ratio = object$ratio
  )
}


#' @export
as.data.frame.spacc_zeta <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    order = x$orders,
    zeta = x$zeta,
    sd = x$zeta_sd,
    ratio = x$ratio
  )
}


#' @export
plot.spacc_zeta <- function(x, type = c("decline", "ratio"),
                             log_scale = FALSE, ...) {
  check_suggests("ggplot2")
  type <- match.arg(type)

  if (type == "decline") {
    df <- data.frame(order = x$orders, zeta = x$zeta, sd = x$zeta_sd)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["order"]], y = .data[["zeta"]])) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[["zeta"]] - .data[["sd"]],
                     ymax = .data[["zeta"]] + .data[["sd"]]),
        fill = "#4CAF50", alpha = 0.2
      ) +
      ggplot2::geom_line(linewidth = 1, color = "#4CAF50") +
      ggplot2::geom_point(size = 2, color = "#4CAF50") +
      ggplot2::labs(
        x = "Order (k)",
        y = expression(zeta[k]),
        title = "Zeta Diversity Decline",
        subtitle = sprintf("%d sites, %s method", x$n_sites, x$method)
      ) +
      ggplot2::theme_minimal(base_size = 12)

    if (log_scale) {
      p <- p + ggplot2::scale_y_log10()
    }

  } else if (type == "ratio") {
    df <- data.frame(order = x$orders, ratio = x$ratio)
    df <- df[!is.na(df$ratio), ]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["order"]], y = .data[["ratio"]])) +
      ggplot2::geom_line(linewidth = 1, color = "#2196F3") +
      ggplot2::geom_point(size = 2, color = "#2196F3") +
      ggplot2::labs(
        x = "Order (k)",
        y = expression(zeta[k] / zeta[k - 1]),
        title = "Zeta Ratio",
        subtitle = "Constant = exponential decline, increasing = power-law"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  }

  p
}
