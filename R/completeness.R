#' Sample Completeness Profile
#'
#' Compute the ratio of observed to estimated diversity across diversity orders,
#' measuring how complete a sample is at each level of the Hill number spectrum.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity to evaluate. Default
#'   `seq(0, 2, by = 0.2)`.
#' @param coords Optional data.frame with columns `x` and `y` for spatial
#'   mapping. When provided, enables `plot(type = "map")` and [as_sf()].
#'
#' @return An object of class `spacc_completeness` containing:
#'   \item{completeness}{Named numeric vector of completeness ratios per q}
#'   \item{observed}{Named numeric vector of observed Hill numbers per q}
#'   \item{estimated}{Named numeric vector of estimated asymptotic Hill numbers per q}
#'   \item{per_site}{Matrix of per-site completeness (sites x q values), or `NULL`}
#'   \item{q}{Vector of diversity orders}
#'   \item{coords}{Coordinates if provided}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Number of species}
#'
#' @details
#' Sample completeness at order q is:
#' \deqn{C_q = \frac{D_q^{obs}}{D_q^{est}}}
#' where \eqn{D_q^{obs}} is the observed Hill number and \eqn{D_q^{est}} is
#' the estimated asymptotic Hill number.
#'
#' Completeness near 1 means the sample captures most of the true diversity
#' at that order. Completeness typically increases with q because dominant
#' species are detected early.
#'
#' Asymptotic estimators used:
#' - q = 0: Chao1 estimator
#' - q = 1: Chao & Jost (2015) entropy estimator
#' - q = 2: Inverse Simpson estimator with bias correction
#' - Other q: Interpolated between adjacent integer estimates
#'
#' When `coords` is provided, per-site completeness is computed by treating
#' each site's abundance vector as an independent sample.
#'
#' @references
#' Chao, A. & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. Ecology, 93,
#' 2533-2547.
#'
#' Chao, A. & Jost, L. (2015). Estimating diversity and entropy profiles
#' via discovery rates of new species. Methods in Ecology and Evolution,
#' 6, 873-882.
#'
#' @seealso [diversityProfile()] for observed diversity profiles,
#'   [chao1()] for richness estimation
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' comp <- completenessProfile(species)
#' print(comp)
#'
#' \donttest{
#' plot(comp)
#' }
#'
#' @export
completenessProfile <- function(x, q = seq(0, 2, by = 0.2),
                                 coords = NULL) {
  x <- as.matrix(x)
  n_sites <- nrow(x)
  n_species <- ncol(x)

  if (!is.null(coords)) {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    stopifnot("x and coords must have same number of rows" = nrow(x) == nrow(coords))
  }

  # Regional (pooled) completeness
  pooled <- colSums(x)
  observed <- sapply(q, function(qi) calc_hill_number(pooled, qi))
  estimated <- sapply(q, function(qi) estimate_asymptotic_hill(pooled, qi))
  completeness <- ifelse(estimated > 0, observed / estimated, 1)
  completeness <- pmin(completeness, 1)  # Cap at 1

  names(observed) <- names(estimated) <- names(completeness) <- paste0("q", q)

  # Per-site completeness
  per_site <- NULL
  if (!is.null(coords)) {
    per_site <- matrix(NA, nrow = n_sites, ncol = length(q))
    colnames(per_site) <- paste0("q", q)
    for (i in seq_len(n_sites)) {
      ab <- as.numeric(x[i, ])
      obs_i <- sapply(q, function(qi) calc_hill_number(ab, qi))
      est_i <- sapply(q, function(qi) estimate_asymptotic_hill(ab, qi))
      per_site[i, ] <- pmin(ifelse(est_i > 0, obs_i / est_i, 1), 1)
    }
  }

  structure(
    list(
      completeness = completeness,
      observed = observed,
      estimated = estimated,
      per_site = per_site,
      q = q,
      coords = coords,
      n_sites = n_sites,
      n_species = n_species
    ),
    class = "spacc_completeness"
  )
}


#' Estimate asymptotic Hill number
#' @noRd
estimate_asymptotic_hill <- function(abundances, q) {
  abundances <- abundances[abundances > 0]
  n <- sum(abundances)
  S_obs <- length(abundances)

  if (S_obs == 0 || n == 0) return(0)

  f1 <- sum(abundances == 1)
  f2 <- sum(abundances == 2)

  if (q == 0) {
    # Chao1 estimator
    if (f2 > 0) {
      return(S_obs + f1^2 / (2 * f2))
    } else {
      return(S_obs + f1 * (f1 - 1) / 2)
    }
  }

  if (abs(q - 1) < 1e-8) {
    # Chao & Jost (2015) entropy estimator
    p <- abundances / n
    # Bias-corrected Shannon entropy
    A <- ifelse(f2 > 0, 2 * f2 / ((n - 1) * f1 + 2 * f2), 0)
    if (f1 == 0 || A == 0) {
      # No correction needed
      H <- -sum(p * log(p))
      return(exp(H))
    }
    # Estimated entropy using coverage-based approach
    H_obs <- -sum(p * log(p))
    # Add contribution from unseen species
    H_unseen <- f1 / n * (1 - A)^(-n + 1) *
      (-log(A) - sum(1 / seq_len(n - 1) * (1 - A)^seq_len(n - 1)))
    H_est <- H_obs + H_unseen
    return(exp(H_est))
  }

  if (abs(q - 2) < 1e-8) {
    # Bias-corrected inverse Simpson
    # D2_est = (n-1) / (sum(a_i*(a_i-1)) / (n*(n-1)))  but Chao-corrected
    sum_p2 <- sum(abundances * (abundances - 1)) / (n * (n - 1))
    if (sum_p2 <= 0) return(S_obs)
    D2_obs <- 1 / sum_p2
    # Chao correction for unseen species
    if (f2 > 0) {
      lambda <- f1 * (f1 - 1) / (2 * (n - 1) * f2)
    } else {
      lambda <- f1 * (f1 - 1) / (2 * (n - 1))
    }
    return(D2_obs + lambda * n / (n - 1))
  }

  # General q: interpolate between bracketing integers
  q_lo <- floor(q)
  q_hi <- ceiling(q)
  if (q_lo == q_hi) q_hi <- q_lo + 1L

  D_lo <- estimate_asymptotic_hill(abundances, q_lo)
  D_hi <- estimate_asymptotic_hill(abundances, q_hi)

  frac <- q - q_lo
  # Log-linear interpolation
  if (D_lo > 0 && D_hi > 0) {
    exp(log(D_lo) * (1 - frac) + log(D_hi) * frac)
  } else {
    D_lo * (1 - frac) + D_hi * frac
  }
}


# S3 methods ---------------------------------------------------------------

#' @export
print.spacc_completeness <- function(x, ...) {
  cat(sprintf("Sample Completeness Profile: %d sites, %d species\n",
              x$n_sites, x$n_species))
  cat(strrep("-", 40), "\n")
  df <- data.frame(
    q = x$q,
    observed = round(x$observed, 2),
    estimated = round(x$estimated, 2),
    completeness = round(x$completeness, 3)
  )
  print(df, row.names = FALSE)
  invisible(x)
}


#' @export
summary.spacc_completeness <- function(object, ...) {
  data.frame(
    q = object$q,
    observed = object$observed,
    estimated = object$estimated,
    completeness = object$completeness
  )
}


#' @export
as.data.frame.spacc_completeness <- function(x, row.names = NULL,
                                              optional = FALSE, ...) {
  data.frame(
    q = x$q,
    observed = x$observed,
    estimated = x$estimated,
    completeness = x$completeness
  )
}


#' @export
plot.spacc_completeness <- function(x, type = c("profile", "comparison", "map"),
                                     point_size = 3, palette = "viridis", ...) {
  check_suggests("ggplot2")
  type <- match.arg(type)

  if (type == "profile") {
    df <- data.frame(q = x$q, completeness = x$completeness)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["q"]], y = .data[["completeness"]])) +
      ggplot2::geom_line(linewidth = 1, color = "#4CAF50") +
      ggplot2::geom_point(size = 2, color = "#4CAF50") +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      ggplot2::scale_y_continuous(limits = c(0, 1.05)) +
      ggplot2::labs(
        x = "Diversity order (q)",
        y = "Completeness (observed / estimated)",
        title = "Sample Completeness Profile",
        subtitle = sprintf("%d sites, %d species", x$n_sites, x$n_species)
      ) +
      ggplot2::theme_minimal(base_size = 12)

  } else if (type == "comparison") {
    df <- data.frame(
      q = rep(x$q, 2),
      value = c(x$observed, x$estimated),
      type = rep(c("Observed", "Estimated"), each = length(x$q))
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["q"]], y = .data[["value"]],
                                           color = .data[["type"]], linetype = .data[["type"]])) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = c("Observed" = "#2196F3", "Estimated" = "#FF9800")) +
      ggplot2::scale_linetype_manual(values = c("Observed" = "solid", "Estimated" = "dashed")) +
      ggplot2::labs(
        x = "Diversity order (q)",
        y = "Effective number of species",
        color = NULL, linetype = NULL,
        title = "Observed vs Estimated Diversity",
        subtitle = sprintf("%d sites, %d species", x$n_sites, x$n_species)
      ) +
      ggplot2::theme_minimal(base_size = 12)

  } else if (type == "map") {
    if (is.null(x$per_site) || is.null(x$coords)) {
      stop("Map requires per-site data. Rerun completenessProfile() with coords.")
    }
    # Map mean completeness across q values
    df <- data.frame(
      x = x$coords$x,
      y = x$coords$y,
      completeness = rowMeans(x$per_site)
    )
    p <- plot_spatial_map(df, "completeness",
                          title = "Mean Sample Completeness",
                          subtitle = sprintf("%d sites, q = [%.1f, %.1f]",
                                             x$n_sites, min(x$q), max(x$q)),
                          point_size = point_size, palette = palette)
  }

  p
}


#' @export
as_sf.spacc_completeness <- function(x, crs = NULL) {
  if (is.null(x$per_site) || is.null(x$coords)) {
    stop("No per-site data. Rerun completenessProfile() with coords.")
  }
  df <- as.data.frame(x$per_site)
  df$completeness_mean <- rowMeans(x$per_site)
  df$x <- x$coords$x
  df$y <- x$coords$y
  as_sf_from_df(df, crs = crs)
}
