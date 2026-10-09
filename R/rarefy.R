#' Individual-Based Rarefaction
#'
#' Compute individual-based rarefaction curves for Hill numbers at any order q.
#' This complements the sample-based accumulation in [spacc()].
#'
#' @param x A site-by-species matrix with abundance data (not presence/absence).
#' @param n_individuals Integer vector. Sample sizes to compute expected
#'   diversity for. Default `NULL` computes for all levels from 1 to total.
#' @param q Numeric. Order of Hill number; any value `>= 0`. Default 0
#'   (species richness). q=1 gives rarefied Shannon diversity, q=2 gives
#'   rarefied Simpson diversity.
#' @param n_boot Integer. Number of bootstrap replicates for CI. Default 100.
#'
#' @return An object of class `spacc_rare` containing:
#'   \item{n}{Sample sizes}
#'   \item{expected}{Expected diversity (Hill number of order q)}
#'   \item{sd}{Standard deviation}
#'   \item{lower, upper}{95 percent confidence bounds}
#'   \item{q}{Order of diversity used}
#'
#' @details
#' The expected Hill number of a random subsample of `n` individuals drawn
#' without replacement from the pooled sample (Chao et al. 2014). Each
#' species' subsample count is hypergeometric, and
#' \deqn{{}^qD(n) = \left(\sum_i E\left[(Y_i / n)^q\right]\right)^{1/(1-q)},}
#' with the exponential of the expected Shannon entropy at q = 1. At q = 0 this
#' is the Hurlbert (1971) formula. Sizes above the sample size are
#' extrapolated towards the asymptotic estimators (Chao et al. 2014; integer
#' `q` only).
#'
#' @examples
#' \donttest{
#' abundance_matrix <- matrix(rpois(50 * 30, 2), nrow = 50)
#' rare <- rarefy(abundance_matrix)
#' plot(rare)
#'
#' # Shannon rarefaction
#' rare_q1 <- rarefy(abundance_matrix, q = 1)
#' plot(rare_q1)
#' }
#'
#' @references
#' Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique
#' and alternative parameters. Ecology, 52, 577-586.
#'
#' Chao, A., Gotelli, N.J., Hsieh, T.C., et al. (2014). Rarefaction and
#' extrapolation with Hill numbers: a framework for sampling and estimation
#' in species diversity studies. Ecological Monographs, 84, 45-67.
#'
#' @export
rarefy <- function(x, n_individuals = NULL, q = 0, n_boot = 100L) {

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

  stopifnot("q must be a single non-negative number" =
              is.numeric(q) && length(q) == 1 && q >= 0)

  rarefy_fn <- function(totals) {
    calc_hill_at_size(as.numeric(totals), as.numeric(n_individuals), q)$hill[1, ]
  }
  expected <- rarefy_fn(species_totals)

  # Bootstrap for variance
  boot_results <- matrix(NA, n_boot, length(n_individuals))
  pool <- rep(seq_along(species_totals), times = species_totals)

  for (b in seq_len(n_boot)) {
    boot_sample <- sample(pool, n_total, replace = TRUE)
    boot_totals <- tabulate(boot_sample, nbins = n_species)
    boot_results[b, ] <- rarefy_fn(boot_totals[boot_totals > 0])
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
      n_boot = n_boot,
      q = q
    ),
    class = "spacc_rare"
  )
}


#' @export
print.spacc_rare <- function(x, ...) {
  q_label <- if (!is.null(x$q) && x$q > 0) sprintf(" (q=%g)", x$q) else ""
  cat(sprintf("Individual-based rarefaction%s\n", q_label))
  cat(strrep("-", 28), "\n")
  cat(sprintf("Total individuals: %d\n", x$n_total))
  cat(sprintf("Observed species: %d\n", x$n_species))
  if (!is.null(x$q) && x$q > 0) cat(sprintf("Diversity order: q=%g\n", x$q))
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

  q_label <- if (!is.null(x$q) && x$q > 0) sprintf(", q=%g", x$q) else ""

  p +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::labs(
      x = "Number of individuals",
      y = if (!is.null(x$q) && x$q > 0) sprintf("Expected diversity (q=%g)", x$q) else "Expected species richness",
      title = "Individual-Based Rarefaction Curve",
      subtitle = sprintf("%d individuals, %d species%s", x$n_total, x$n_species, q_label)
    ) +
    spacc_theme()
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
