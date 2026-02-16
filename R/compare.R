#' Compare Two Accumulation Curves
#'
#' Test whether two species accumulation curves differ significantly.
#'
#' @param x A `spacc` object.
#' @param y A `spacc` object.
#' @param method Character. Comparison method: `"permutation"` (default),
#'   `"bootstrap"`, or `"auc"` (area under curve difference).
#' @param normalize Logical. If `TRUE`, each seed's curve is divided by its
#'   final value before computing AUC, so that curves are compared on a
#'   \[0, 1\] scale (proportion of total species). This compares the **shape**
#'   of accumulation rather than absolute species counts. Default `FALSE`.
#' @param n_perm Integer. Number of permutations/bootstrap replicates. Default 999.
#' @param ... Additional arguments passed to comparison methods.
#'
#' @return An object of class `spacc_comp` containing:
#'   \item{x_name, y_name}{Names of compared objects}
#'   \item{auc_diff}{Difference in area under curve}
#'   \item{p_value}{P-value from permutation test}
#'   \item{saturation_diff}{Difference in saturation points}
#'   \item{method}{Comparison method used}
#'   \item{normalized}{Whether curves were normalized before comparison}
#'
#' @examples
#' \dontrun{
#' sac_native <- spacc(native_species, coords)
#' sac_alien <- spacc(alien_species, coords)
#'
#' # Compare absolute accumulation
#' comp <- compare(sac_native, sac_alien)
#'
#' # Compare curve shapes (normalized to [0, 1])
#' comp_norm <- compare(sac_native, sac_alien, normalize = TRUE)
#' print(comp_norm)
#' plot(comp_norm)
#' }
#'
#' @references
#' Colwell, R.K., Mao, C.X. & Chang, J. (2004). Interpolating, extrapolating,
#' and comparing incidence-based species accumulation curves. Ecology, 85,
#' 2717-2727.
#'
#' Gotelli, N.J. & Colwell, R.K. (2011). Estimating species richness. In:
#' Biological Diversity: Frontiers in Measurement and Assessment, pp. 39-54.
#' Oxford University Press.
#'
#' @export
compare <- function(x, y, method = c("permutation", "bootstrap", "auc"),
                    normalize = FALSE, n_perm = 999L, ...) {

  method <- match.arg(method)

  stopifnot(
    "x must be a spacc object" = inherits(x, "spacc"),
    "y must be a spacc object" = inherits(y, "spacc"),
    "x and y must have same number of sites" = x$n_sites == y$n_sites
  )

  # Get names from call
  x_name <- deparse(substitute(x))
  y_name <- deparse(substitute(y))

  # Normalize: divide each seed's curve by its final value -> [0, 1]
  curves_x <- x$curves
  curves_y <- y$curves
  if (normalize) {
    row_max_x <- curves_x[, ncol(curves_x)]
    row_max_y <- curves_y[, ncol(curves_y)]
    row_max_x[row_max_x == 0] <- 1
    row_max_y[row_max_y == 0] <- 1
    curves_x <- curves_x / row_max_x
    curves_y <- curves_y / row_max_y
  }

  # Area under curve (sum of mean curve)
  auc_x <- sum(colMeans(curves_x))
  auc_y <- sum(colMeans(curves_y))
  auc_diff <- auc_x - auc_y

  # Compute summary stats for saturation
  summ_x <- summary(x)
  summ_y <- summary(y)
  sat_diff <- summ_x$saturation_point - summ_y$saturation_point

  # Permutation test for AUC difference
  if (method == "permutation") {
    all_curves <- rbind(curves_x, curves_y)
    n_x <- nrow(curves_x)
    n_total <- nrow(all_curves)

    null_diffs <- vapply(seq_len(n_perm), function(i) {
      perm_idx <- sample(n_total)
      perm_x <- all_curves[perm_idx[1:n_x], , drop = FALSE]
      perm_y <- all_curves[perm_idx[(n_x + 1):n_total], , drop = FALSE]
      sum(colMeans(perm_x)) - sum(colMeans(perm_y))
    }, numeric(1))

    p_value <- mean(abs(null_diffs) >= abs(auc_diff))
  } else if (method == "bootstrap") {
    boot_x <- vapply(seq_len(n_perm), function(i) {
      idx <- sample(nrow(curves_x), replace = TRUE)
      sum(colMeans(curves_x[idx, , drop = FALSE]))
    }, numeric(1))

    boot_y <- vapply(seq_len(n_perm), function(i) {
      idx <- sample(nrow(curves_y), replace = TRUE)
      sum(colMeans(curves_y[idx, , drop = FALSE]))
    }, numeric(1))

    boot_diff <- boot_x - boot_y
    ci <- stats::quantile(boot_diff, c(0.025, 0.975))
    p_value <- ifelse(ci[1] > 0 || ci[2] < 0, 0.05, 0.5) # rough approximation
  } else {
    # Simple AUC comparison, no p-value
    p_value <- NA
  }

  structure(
    list(
      x_name = x_name,
      y_name = y_name,
      x = x,
      y = y,
      auc_x = auc_x,
      auc_y = auc_y,
      auc_diff = auc_diff,
      saturation_x = summ_x$saturation_point,
      saturation_y = summ_y$saturation_point,
      saturation_diff = sat_diff,
      p_value = p_value,
      n_perm = n_perm,
      method = method,
      normalized = normalize
    ),
    class = "spacc_comp"
  )
}


#' @export
print.spacc_comp <- function(x, ...) {
  cat("Comparison:", x$x_name, "vs", x$y_name, "\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("Method: %s (n=%d)\n", x$method, x$n_perm))
  if (isTRUE(x$normalized)) cat("Normalized: yes (shape comparison)\n")
  cat(sprintf("AUC difference: %.1f", x$auc_diff))
  if (!is.na(x$p_value)) {
    sig <- ifelse(x$p_value < 0.001, "***",
           ifelse(x$p_value < 0.01, "**",
           ifelse(x$p_value < 0.05, "*", "")))
    p_str <- if (x$p_value < 0.001) "p < 0.001" else sprintf("p = %.3f", x$p_value)
    cat(sprintf(" (%s%s)\n", p_str, sig))
  } else {
    cat("\n")
  }
  cat(sprintf("Saturation: %s at %d sites, %s at %d sites\n",
              x$x_name, x$saturation_x,
              x$y_name, x$saturation_y))

  # Interpretation
  if (!is.na(x$p_value) && x$p_value < 0.05) {
    if (x$saturation_x < x$saturation_y) {
      cat(sprintf("\n%s saturates faster.\n", x$x_name))
    } else if (x$saturation_x > x$saturation_y) {
      cat(sprintf("\n%s saturates faster.\n", x$y_name))
    }
  }
  invisible(x)
}


#' @export
summary.spacc_comp <- function(object, ...) {
  print(object)
  invisible(object)
}


#' @export
plot.spacc_comp <- function(x, ci = TRUE, ci_alpha = 0.2, ...) {
  check_suggests("ggplot2")

  grouped <- c(x$x, x$y)
  grouped$group_names <- c(x$x_name, x$y_name)
  names(grouped$curves) <- grouped$group_names
  names(grouped$n_species) <- grouped$group_names

  p <- plot(grouped, ci = ci, ci_alpha = ci_alpha, ...)

  # Add significance annotation
  if (!is.na(x$p_value)) {
    sig_text <- if (x$p_value < 0.001) "p < 0.001" else sprintf("p = %.3f", x$p_value)
    p <- p + ggplot2::annotate(
      "text", x = x$x$n_sites * 0.8, y = max(x$x$n_species, x$y$n_species) * 0.1,
      label = sig_text, size = 4
    )
  }

  p + ggplot2::labs(title = sprintf("Comparison: %s vs %s", x$x_name, x$y_name))
}
