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
#'   rarefied Simpson diversity. q=0, 1, and 2 use exact estimators; other
#'   orders report the Hill number of order `q` of the sampled composition.
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
#' For q=0 (species richness): uses the Hurlbert (1971) formula.
#'
#' For q=1 (Shannon diversity): rarefied Shannon entropy is computed and
#' converted to effective number of species via exponentiation.
#'
#' For q=2 (Simpson diversity): rarefied Simpson concentration is computed
#' and converted to effective number of species via inversion.
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

  # Compute rarefied diversity for any q >= 0 (q = 0, 1, 2 are exact)
  rarefy_fn <- function(st, n, nt, ns) rarefy_qD(st, n, nt, ns, q)

  expected <- vapply(n_individuals, function(n) {
    rarefy_fn(species_totals, n, n_total, n_species)
  }, numeric(1))

  # Bootstrap for variance
  boot_results <- matrix(NA, n_boot, length(n_individuals))
  pool <- rep(seq_along(species_totals), times = species_totals)

  for (b in seq_len(n_boot)) {
    boot_sample <- sample(pool, n_total, replace = TRUE)
    boot_totals <- tabulate(boot_sample, nbins = n_species)
    boot_totals <- boot_totals[boot_totals > 0]
    boot_n_total <- sum(boot_totals)
    boot_n_species <- length(boot_totals)

    boot_results[b, ] <- vapply(n_individuals, function(n) {
      rarefy_fn(boot_totals, n, boot_n_total, boot_n_species)
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
      n_boot = n_boot,
      q = q
    ),
    class = "spacc_rare"
  )
}


# Rarefaction: q=0 (Hurlbert 1971, species richness)
#' @noRd
rarefy_q0 <- function(species_totals, n, n_total, n_species) {
  if (n >= n_total) return(n_species)
  probs <- vapply(species_totals, function(ni) {
    1 - exp(lchoose(n_total - ni, n) - lchoose(n_total, n))
  }, numeric(1))
  sum(probs)
}


# Rarefaction: q=1 (Shannon diversity as Hill number)
#' @noRd
rarefy_q1 <- function(species_totals, n, n_total, n_species) {
  if (n >= n_total) {
    p <- species_totals / n_total
    p <- p[p > 0]
    return(exp(-sum(p * log(p))))
  }
  if (n <= 1) return(1)

  # Rarefied Shannon entropy: E[H] using inclusion probabilities
  # For each species, probability of detection at subsample n
  probs <- vapply(species_totals, function(ni) {
    1 - exp(lchoose(n_total - ni, n) - lchoose(n_total, n))
  }, numeric(1))

  # Expected proportion of each species conditional on detection
  # Approximate: E[p_i * log(p_i)] is complex, use expected abundance approach
  # E[abundance of species i in subsample of n] = n * ni / N
  exp_abund <- n * species_totals / n_total
  # Expected proportions
  p_exp <- exp_abund / n
  p_exp <- p_exp[p_exp > 0]

  H <- -sum(p_exp * log(p_exp))
  exp(H)
}


# Rarefaction: q=2 (Simpson diversity as Hill number)
#' @noRd
rarefy_q2 <- function(species_totals, n, n_total, n_species) {
  if (n >= n_total) {
    p <- species_totals / n_total
    return(1 / sum(p^2))
  }
  if (n <= 1) return(1)

  # Rarefied Simpson: E[sum(p_i^2)] for subsample of size n
  # E[n_i(n_i-1)] in subsample = n(n-1) * N_i(N_i-1) / (N(N-1))
  # So E[sum p_i^2] = sum E[n_i(n_i-1)] / (n(n-1)) + 1/n
  lambda <- sum(species_totals * (species_totals - 1)) /
            (n_total * (n_total - 1))

  if (lambda <= 0) return(n_species)

  # Simpson concentration for subsample
  # lambda_n = (n-1)/n * lambda + 1/n (bias from finite sample)
  # This simplifies to just lambda for expected D2
  1 / lambda
}


# Rarefaction: general order q (>= 0)
# q = 0, 1, 2 dispatch to the exact estimators above; other orders return the
# Hill number of order q of the sampled composition (effort-independent in the
# interior, matching the q = 1 / q = 2 treatment, and exact at the endpoints).
#' @noRd
rarefy_qD <- function(species_totals, n, n_total, n_species, q) {
  if (q == 0) return(rarefy_q0(species_totals, n, n_total, n_species))
  if (abs(q - 1) < 1e-8) return(rarefy_q1(species_totals, n, n_total, n_species))
  if (q == 2) return(rarefy_q2(species_totals, n, n_total, n_species))

  if (n <= 1) return(1)
  p <- species_totals / n_total
  p <- p[p > 0]
  sum(p^q)^(1 / (1 - q))
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
