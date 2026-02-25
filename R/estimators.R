#' Chao1 Richness Estimator
#'
#' Estimate total species richness from abundance data using the Chao1
#' estimator (Chao 1984). Uses the number of singletons and doubletons
#' to estimate undetected species.
#'
#' @param x A site-by-species matrix (abundance data). Columns are pooled
#'   across sites.
#'
#' @return An object of class `spacc_estimate` with components:
#' \describe{
#'   \item{estimator}{Name of the estimator (`"chao1"`)}
#'   \item{estimate}{Estimated total richness}
#'   \item{se}{Standard error of the estimate}
#'   \item{lower}{Lower 95 percent confidence bound}
#'   \item{upper}{Upper 95 percent confidence bound}
#'   \item{S_obs}{Observed species richness}
#'   \item{details}{List with `f1` (singletons) and `f2` (doubletons)}
#' }
#'
#' @details
#' The Chao1 estimator is:
#' \deqn{S_{Chao1} = S_{obs} + \frac{f_1^2}{2 f_2}}
#' where \eqn{f_1} is the number of singletons (species with total abundance 1)
#' and \eqn{f_2} is the number of doubletons (abundance 2).
#'
#' When \eqn{f_2 = 0}, the bias-corrected form is used:
#' \deqn{S_{Chao1} = S_{obs} + \frac{f_1 (f_1 - 1)}{2}}
#'
#' @references
#' Chao, A. (1984). Nonparametric estimation of the number of classes in a
#' population. Scandinavian Journal of Statistics, 11, 265-270.
#'
#' Chao, A. (1987). Estimating the population size for capture-recapture data
#' with unequal catchability. Biometrics, 43, 783-791.
#'
#' @seealso [chao2()] for incidence-based estimation, [ace()] for
#'   abundance-based coverage estimator
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' chao1(species)
#'
#' @export
chao1 <- function(x) {
  x <- as.matrix(x)
  abundances <- colSums(x)
  abundances <- abundances[abundances > 0]

  S_obs <- length(abundances)
  f1 <- sum(abundances == 1)
  f2 <- sum(abundances == 2)

  if (f2 > 0) {
    estimate <- S_obs + f1^2 / (2 * f2)
    # Variance (Chao 1987)
    se <- sqrt(f2 * ((f1 / f2)^2 / 4 + (f1 / f2)^3 + (f1 / f2)^4 / 4))
  } else {
    # Bias-corrected form when f2 = 0
    estimate <- S_obs + f1 * (f1 - 1) / 2
    se <- sqrt(f1 * (f1 - 1) / 2 + f1 * (2 * f1 - 1)^2 / 4 - f1^4 / (4 * estimate))
  }

  # Log-transformed CI (Chao 1987)
  ci <- chao_log_ci(estimate, se, S_obs)

  structure(
    list(
      estimator = "chao1",
      estimate = estimate,
      se = se,
      lower = ci[1],
      upper = ci[2],
      S_obs = S_obs,
      details = list(f1 = f1, f2 = f2)
    ),
    class = "spacc_estimate"
  )
}


#' Chao2 Richness Estimator
#'
#' Estimate total species richness from incidence (presence/absence) data
#' using the Chao2 estimator (Chao 1987).
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#'   Automatically binarized.
#'
#' @return An object of class `spacc_estimate`.
#'
#' @details
#' The Chao2 estimator is the incidence-based analogue of Chao1:
#' \deqn{S_{Chao2} = S_{obs} + \frac{Q_1^2}{2 Q_2}}
#' where \eqn{Q_1} is the number of uniques (species found at exactly 1 site)
#' and \eqn{Q_2} is the number of duplicates (species found at exactly 2 sites).
#'
#' @references
#' Chao, A. (1987). Estimating the population size for capture-recapture data
#' with unequal catchability. Biometrics, 43, 783-791.
#'
#' @seealso [chao1()] for abundance-based estimation
#'
#' @examples
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' chao2(species)
#'
#' @export
chao2 <- function(x) {
  x <- as.matrix(x)
  pa <- (x > 0) * 1L

  freq <- colSums(pa)
  freq <- freq[freq > 0]

  S_obs <- length(freq)
  Q1 <- sum(freq == 1)
  Q2 <- sum(freq == 2)
  n <- nrow(pa)

  if (Q2 > 0) {
    estimate <- S_obs + (n - 1) / n * Q1^2 / (2 * Q2)
    se <- sqrt(Q2 * (0.5 * ((n - 1) / n) * (Q1 / Q2)^2 +
                      ((n - 1) / n)^2 * (Q1 / Q2)^3 +
                      0.25 * ((n - 1) / n)^2 * (Q1 / Q2)^4))
  } else {
    estimate <- S_obs + (n - 1) / n * Q1 * (Q1 - 1) / 2
    se <- sqrt((n - 1) / n * Q1 * (Q1 - 1) / 2 +
               ((n - 1) / n)^2 * Q1 * (2 * Q1 - 1)^2 / 4 -
               ((n - 1) / n)^2 * Q1^4 / (4 * max(estimate, 1)))
  }

  ci <- chao_log_ci(estimate, se, S_obs)

  structure(
    list(
      estimator = "chao2",
      estimate = estimate,
      se = se,
      lower = ci[1],
      upper = ci[2],
      S_obs = S_obs,
      details = list(Q1 = Q1, Q2 = Q2, n_sites = n)
    ),
    class = "spacc_estimate"
  )
}


#' ACE Richness Estimator
#'
#' Abundance-based Coverage Estimator (ACE) of total species richness
#' (Chao & Lee 1992). Separates species into rare and abundant groups
#' based on an abundance threshold.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param threshold Integer. Abundance threshold separating rare from abundant
#'   species. Default 10.
#'
#' @return An object of class `spacc_estimate`.
#'
#' @details
#' ACE partitions species into rare (abundance \eqn{\le} threshold) and
#' abundant (abundance > threshold) groups. The estimate is:
#' \deqn{S_{ACE} = S_{abund} + \frac{S_{rare}}{C_{ACE}} + \frac{f_1}{C_{ACE}} \gamma^2}
#' where \eqn{C_{ACE}} is the sample coverage of rare species and
#' \eqn{\gamma^2} is the estimated coefficient of variation.
#'
#' @references
#' Chao, A. & Lee, S.M. (1992). Estimating the number of classes via sample
#' coverage. Journal of the American Statistical Association, 87, 210-217.
#'
#' @seealso [chao1()], [chao2()]
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' ace(species)
#'
#' @export
ace <- function(x, threshold = 10L) {
  x <- as.matrix(x)
  abundances <- colSums(x)
  abundances <- abundances[abundances > 0]

  S_obs <- length(abundances)

  # Separate rare and abundant species
  rare <- abundances[abundances <= threshold]
  S_rare <- length(rare)
  S_abund <- sum(abundances > threshold)

  if (S_rare == 0) {
    # No rare species: ACE = S_obs
    return(structure(
      list(
        estimator = "ace",
        estimate = S_obs,
        se = 0,
        lower = S_obs,
        upper = S_obs,
        S_obs = S_obs,
        details = list(S_rare = 0, S_abund = S_abund,
                       C_ace = 1, gamma_sq = 0, threshold = threshold)
      ),
      class = "spacc_estimate"
    ))
  }

  # Total individuals in rare group
  n_rare <- sum(rare)
  f1 <- sum(abundances == 1)

  # Sample coverage of rare species
  C_ace <- 1 - f1 / n_rare

  if (C_ace <= 0) {
    # Degenerate case: all rare species are singletons
    estimate <- S_obs + f1 * (f1 - 1) / 2
    se <- sqrt(f1 * (f1 - 1) / 2)
    ci <- chao_log_ci(estimate, se, S_obs)

    return(structure(
      list(
        estimator = "ace",
        estimate = estimate,
        se = se,
        lower = ci[1],
        upper = ci[2],
        S_obs = S_obs,
        details = list(S_rare = S_rare, S_abund = S_abund,
                       C_ace = C_ace, gamma_sq = NA, threshold = threshold)
      ),
      class = "spacc_estimate"
    ))
  }

  # Coefficient of variation squared
  freq_counts <- tabulate(rare, nbins = threshold)
  i_vals <- seq_len(threshold)
  gamma_sq <- max(0, S_rare / C_ace *
    sum(i_vals * (i_vals - 1) * freq_counts) / (n_rare * (n_rare - 1)) - 1)

  estimate <- S_abund + S_rare / C_ace + f1 / C_ace * gamma_sq

  # Approximate SE
  se <- sqrt(sum((seq_len(threshold))^2 * freq_counts *
                  (1 - freq_counts / S_rare)) / (C_ace^2))

  ci <- chao_log_ci(estimate, se, S_obs)

  structure(
    list(
      estimator = "ace",
      estimate = estimate,
      se = se,
      lower = ci[1],
      upper = ci[2],
      S_obs = S_obs,
      details = list(S_rare = S_rare, S_abund = S_abund,
                     C_ace = C_ace, gamma_sq = gamma_sq, threshold = threshold)
    ),
    class = "spacc_estimate"
  )
}


#' Jackknife Richness Estimator
#'
#' First- or second-order jackknife estimator of total species richness
#' (Burnham & Overton 1978, 1979).
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#'   Automatically binarized.
#' @param order Integer. Jackknife order: 1 (default) or 2.
#'
#' @return An object of class `spacc_estimate`.
#'
#' @details
#' First-order jackknife:
#' \deqn{S_{jack1} = S_{obs} + Q_1 \frac{n-1}{n}}
#'
#' Second-order jackknife:
#' \deqn{S_{jack2} = S_{obs} + \frac{Q_1(2n-3)}{n} - \frac{Q_2(n-2)^2}{n(n-1)}}
#'
#' @references
#' Burnham, K.P. & Overton, W.S. (1978). Estimation of the size of a closed
#' population when capture probabilities vary among animals. Biometrika, 65,
#' 625-633.
#'
#' Burnham, K.P. & Overton, W.S. (1979). Robust estimation of population size
#' when capture probabilities vary among animals. Ecology, 60, 927-936.
#'
#' @seealso [chao2()], [bootstrap_richness()]
#'
#' @examples
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' jackknife(species, order = 1)
#' jackknife(species, order = 2)
#'
#' @export
jackknife <- function(x, order = 1L) {
  x <- as.matrix(x)
  pa <- (x > 0) * 1L

  if (!order %in% c(1L, 2L)) {
    stop("order must be 1 or 2")
  }

  freq <- colSums(pa)
  freq <- freq[freq > 0]

  S_obs <- length(freq)
  n <- nrow(pa)
  Q1 <- sum(freq == 1)
  Q2 <- sum(freq == 2)

  if (order == 1L) {
    estimate <- S_obs + Q1 * (n - 1) / n
    # Variance (Smith & van Belle 1984)
    se <- sqrt(sum((n - freq)^2 * (freq == 1)) * (n - 1) / n -
               Q1^2 / n)
    # Simplified SE if above fails
    if (is.na(se) || se < 0) se <- sqrt(Q1 * (n - 1)^2 / n^2)
  } else {
    estimate <- S_obs + Q1 * (2 * n - 3) / n - Q2 * (n - 2)^2 / (n * (n - 1))
    se <- sqrt(Q1 * (2 * n - 3)^2 / n^2 +
               Q2 * (n - 2)^4 / (n * (n - 1))^2 -
               (Q1 * (2 * n - 3) / n - Q2 * (n - 2)^2 / (n * (n - 1)))^2 / S_obs)
    if (is.na(se) || se < 0) se <- sqrt(Q1)
  }

  ci <- c(estimate - 1.96 * se, estimate + 1.96 * se)
  ci[1] <- max(S_obs, ci[1])

  structure(
    list(
      estimator = paste0("jackknife", order),
      estimate = estimate,
      se = se,
      lower = ci[1],
      upper = ci[2],
      S_obs = S_obs,
      details = list(Q1 = Q1, Q2 = Q2, n_sites = n, order = order)
    ),
    class = "spacc_estimate"
  )
}


#' Bootstrap Richness Estimator
#'
#' Bootstrap estimator of total species richness (Smith & van Belle 1984).
#' Uses species detection probabilities to estimate undetected species.
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#'   Automatically binarized.
#' @param n_boot Integer. Number of bootstrap replicates for SE. Default 200.
#'
#' @return An object of class `spacc_estimate`.
#'
#' @details
#' The bootstrap estimator is:
#' \deqn{S_{boot} = S_{obs} + \sum_{i=1}^{S_{obs}} (1 - p_i)^n}
#' where \eqn{p_i} is the proportion of sites where species \eqn{i} occurs.
#'
#' @references
#' Smith, E.P. & van Belle, G. (1984). Nonparametric estimation of species
#' richness. Biometrics, 40, 119-129.
#'
#' @seealso [jackknife()], [chao2()]
#'
#' @examples
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' bootstrap_richness(species, n_boot = 100)
#'
#' @export
bootstrap_richness <- function(x, n_boot = 200L) {
  x <- as.matrix(x)
  pa <- (x > 0) * 1L

  n <- nrow(pa)
  freq <- colSums(pa)
  freq <- freq[freq > 0]

  S_obs <- length(freq)
  p <- freq / n

  # Point estimate
  estimate <- S_obs + sum((1 - p)^n)

  # Bootstrap SE
  boot_estimates <- vapply(seq_len(n_boot), function(b) {
    idx <- sample(n, replace = TRUE)
    boot_pa <- pa[idx, , drop = FALSE]
    boot_freq <- colSums(boot_pa)
    boot_freq <- boot_freq[boot_freq > 0]
    boot_S <- length(boot_freq)
    boot_p <- boot_freq / n
    boot_S + sum((1 - boot_p)^n)
  }, numeric(1))

  se <- stats::sd(boot_estimates)
  ci <- stats::quantile(boot_estimates, c(0.025, 0.975))

  structure(
    list(
      estimator = "bootstrap",
      estimate = estimate,
      se = se,
      lower = max(S_obs, unname(ci[1])),
      upper = unname(ci[2]),
      S_obs = S_obs,
      details = list(n_boot = n_boot, n_sites = n)
    ),
    class = "spacc_estimate"
  )
}


# Log-transformed CI for Chao estimators (Chao 1987)
chao_log_ci <- function(estimate, se, S_obs) {
  T_val <- estimate - S_obs
  if (T_val <= 0 || se <= 0) return(c(S_obs, estimate))

  K <- exp(1.96 * sqrt(log(1 + (se / T_val)^2)))
  c(S_obs + T_val / K, S_obs + T_val * K)
}


# S3 methods for spacc_estimate -----------------------------------------

#' @export
print.spacc_estimate <- function(x, ...) {
  cat(sprintf("Richness Estimator: %s\n", x$estimator))
  cat(strrep("-", 32), "\n")
  cat(sprintf("Observed species (S_obs): %d\n", x$S_obs))
  cat(sprintf("Estimated richness:       %.1f\n", x$estimate))
  cat(sprintf("Standard error:           %.1f\n", x$se))
  cat(sprintf("95%% CI:                   [%.1f, %.1f]\n", x$lower, x$upper))
  invisible(x)
}


#' @export
summary.spacc_estimate <- function(object, ...) {
  df <- data.frame(
    estimator = object$estimator,
    S_obs = object$S_obs,
    estimate = object$estimate,
    se = object$se,
    lower = object$lower,
    upper = object$upper
  )
  df
}


#' @export
as.data.frame.spacc_estimate <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    estimator = x$estimator,
    S_obs = x$S_obs,
    estimate = x$estimate,
    se = x$se,
    lower = x$lower,
    upper = x$upper
  )
}


#' @export
plot.spacc_estimate <- function(x, ...) {
  check_suggests("ggplot2")

  df <- data.frame(
    label = c("Observed", x$estimator),
    value = c(x$S_obs, x$estimate),
    lower = c(x$S_obs, x$lower),
    upper = c(x$S_obs, x$upper)
  )
  df$label <- factor(df$label, levels = df$label)

  ggplot2::ggplot(df, ggplot2::aes(x = .data[["label"]], y = .data[["value"]])) +
    ggplot2::geom_col(fill = c("#78909C", "#4CAF50"), width = 0.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]),
      width = 0.15, linewidth = 0.8
    ) +
    ggplot2::labs(
      x = NULL, y = "Species richness",
      title = sprintf("Richness Estimation: %s", x$estimator)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
