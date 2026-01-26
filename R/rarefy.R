#' Individual-Based Rarefaction
#'
#' Compute classic individual-based rarefaction curves. This complements
#' the sample-based accumulation in [spacc()].
#'
#' @param x A site-by-species matrix with abundance data (not presence/absence).
#' @param n_individuals Integer vector. Sample sizes to compute expected
#'   richness for. Default `NULL` computes for all levels from 1 to total.
#' @param n_boot Integer. Number of bootstrap replicates for CI. Default 100.
#'
#' @return An object of class `spacc_rare` containing:
#'   \item{n}{Sample sizes}
#'   \item{expected}{Expected species richness}
#'   \item{sd}{Standard deviation}
#'   \item{lower, upper}{95% confidence bounds}
#'
#' @examples
#' \dontrun{
#' # With abundance data
#' rare <- rarefy(abundance_matrix)
#' plot(rare)
#' }
#'
#' @export
rarefy <- function(x, n_individuals = NULL, n_boot = 100L) {

  x <- as.matrix(x)

  # Total individuals per species
  species_totals <- colSums(x)
  species_totals <- species_totals[species_totals > 0]

  n_total <- sum(species_totals)
  n_species <- length(species_totals)

  # Sample sizes to compute
 if (is.null(n_individuals)) {
    n_individuals <- seq(1, n_total, length.out = min(100, n_total))
    n_individuals <- unique(round(n_individuals))
  }

  # Analytical expected richness (Hurlbert 1971)
  expected <- vapply(n_individuals, function(n) {
    if (n >= n_total) return(n_species)
    # E[S] = sum(1 - choose(N-Ni, n) / choose(N, n))
    probs <- vapply(species_totals, function(ni) {
      1 - exp(lchoose(n_total - ni, n) - lchoose(n_total, n))
    }, numeric(1))
    sum(probs)
  }, numeric(1))

  # Bootstrap for variance
  boot_results <- matrix(NA, n_boot, length(n_individuals))

  # Create pool of individuals
  pool <- rep(seq_along(species_totals), times = species_totals)

  for (b in seq_len(n_boot)) {
    boot_sample <- sample(pool, n_total, replace = TRUE)
    boot_totals <- tabulate(boot_sample, nbins = n_species)
    boot_totals <- boot_totals[boot_totals > 0]
    boot_n_total <- sum(boot_totals)

    boot_results[b, ] <- vapply(n_individuals, function(n) {
      if (n >= boot_n_total) return(length(boot_totals))
      probs <- vapply(boot_totals, function(ni) {
        1 - exp(lchoose(boot_n_total - ni, n) - lchoose(boot_n_total, n))
      }, numeric(1))
      sum(probs)
    }, numeric(1))
  }

  sd_vals <- apply(boot_results, 2, stats::sd, na.rm = TRUE)
  lower <- apply(boot_results, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(boot_results, 2, stats::quantile, probs = 0.975, na.rm = TRUE)

  structure(
    list(
      n = n_individuals,
      expected = expected,
      sd = sd_vals,
      lower = lower,
      upper = upper,
      n_total = n_total,
      n_species = n_species,
      n_boot = n_boot
    ),
    class = "spacc_rare"
  )
}


#' @export
print.spacc_rare <- function(x, ...) {
  cat("Individual-based rarefaction\n")
  cat(strrep("-", 28), "\n")
  cat(sprintf("Total individuals: %d\n", x$n_total))
  cat(sprintf("Observed species: %d\n", x$n_species))
  cat(sprintf("Bootstrap replicates: %d\n", x$n_boot))
  invisible(x)
}


#' @export
plot.spacc_rare <- function(x, ci = TRUE, ci_alpha = 0.3, ...) {
  check_suggests("ggplot2")

  df <- data.frame(
    n = x$n,
    expected = x$expected,
    lower = x$lower,
    upper = x$upper
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = n, y = expected))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::labs(
      x = "Number of individuals",
      y = "Expected species richness",
      title = "Individual-Based Rarefaction Curve",
      subtitle = sprintf("%d individuals, %d species", x$n_total, x$n_species)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' @export
as.data.frame.spacc_rare <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    n = x$n,
    expected = x$expected,
    sd = x$sd,
    lower = x$lower,
    upper = x$upper
  )
}
