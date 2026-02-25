#' Alpha Diversity (Per-Site)
#'
#' Compute Hill numbers for each site individually.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#' @param coords Optional data.frame with columns `x` and `y` for spatial mapping.
#'   When provided, returns a `spacc_alpha` object with [as_sf()] and
#'   `plot(type = "map")` support.
#'
#' @return If `coords` is `NULL`, a data.frame with columns for each q value.
#'   If `coords` is provided, a `spacc_alpha` object.
#'
#' @details
#' Alpha diversity represents local (within-site) diversity. For Hill numbers:
#' - q=0: Species richness
#' - q=1: Exponential of Shannon entropy
#' - q=2: Inverse Simpson concentration
#'
#' @references
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta
#' components. Ecology, 88, 2427-2439.
#'
#' @seealso [gammaDiversity()] for regional diversity, [diversityPartition()]
#'   for full alpha-beta-gamma decomposition
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' alpha <- alphaDiversity(species, q = c(0, 1, 2))
#' head(alpha)
#'
#' # Mean alpha diversity
#' colMeans(alpha)
#'
#' @export
alphaDiversity <- function(x, q = c(0, 1, 2), coords = NULL) {
  x <- as.matrix(x)
  n_sites <- nrow(x)
  result <- matrix(NA, nrow = n_sites, ncol = length(q))
  colnames(result) <- paste0("q", q)
  for (i in seq_len(n_sites)) {
    abundances <- as.numeric(x[i, ])
    for (j in seq_along(q)) {
      result[i, j] <- calc_hill_number(abundances, q[j])
    }
  }
  result_df <- as.data.frame(result)

  if (is.null(coords)) return(result_df)

  stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
  stopifnot("x and coords must have same number of rows" = nrow(x) == nrow(coords))

  result_df$x <- coords$x
  result_df$y <- coords$y

  structure(
    list(
      values = result_df,
      q = q,
      coords = coords,
      n_sites = n_sites,
      n_species = ncol(x)
    ),
    class = "spacc_alpha"
  )
}
#' Gamma Diversity (Regional)
#'
#' Compute Hill numbers for the pooled community across all sites.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#'
#' @return A named numeric vector with gamma diversity for each q.
#'
#' @details
#' Gamma diversity represents regional (total) diversity across all sites.
#' It is computed by pooling abundances across all sites and calculating
#' Hill numbers on the combined community.
#'
#' @references
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta
#' components. Ecology, 88, 2427-2439.
#'
#' @seealso [alphaDiversity()] for local diversity, [diversityPartition()]
#'   for full alpha-beta-gamma decomposition
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' gammaDiversity(species, q = c(0, 1, 2))
#'
#' @export
gammaDiversity <- function(x, q = c(0, 1, 2)) {
  x <- as.matrix(x)
  # Pool abundances across all sites
  pooled <- colSums(x)
  result <- sapply(q, function(qi) calc_hill_number(pooled, qi))
  names(result) <- paste0("q", q)
  result
}
#' Alpha-Beta-Gamma Diversity Partitioning
#'
#' Decompose regional (gamma) diversity into local (alpha) and turnover (beta)
#' components using multiplicative partitioning of Hill numbers.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#' @param weights Character or numeric. Site weights for alpha calculation:
#'   - `"equal"`: Equal weights (default)
#'   - `"proportional"`: Weights proportional to site abundance
#'   - Numeric vector of custom weights
#' @param coords Optional data.frame with columns `x` and `y` for spatial mapping.
#'   When provided, enables [as_sf()] and `plot(type = "map")`.
#'
#' @return An object of class `spacc_partition` containing:
#'   \item{alpha}{Mean alpha diversity (effective number of species per site)}
#'   \item{beta}{Beta diversity (effective number of communities)}
#'   \item{gamma}{Gamma diversity (regional species pool)}
#'   \item{q}{Orders of diversity}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Total species count}
#'
#' @details
#' This function implements Jost (2007) multiplicative partitioning:
#'
#' \deqn{\gamma = \alpha \times \beta}
#'
#' Where:
#' - **Alpha**: Mean effective number of species per site
#' - **Beta**: Effective number of distinct communities (1 = all identical,
#'   n_sites = all completely different)
#' - **Gamma**: Total effective number of species in the region
#'
#' Beta diversity is interpreted as the effective number of communities:
#' - Beta = 1: All sites have identical composition
#' - Beta = n_sites: Sites share no species
#'
#' @references
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta
#' components. Ecology, 88, 2427-2439.
#'
#' Chao, A., Chiu, C.H. & Jost, L. (2014). Unifying species diversity, phylogenetic
#' diversity, functional diversity, and related similarity and differentiation
#' measures through Hill numbers. Annual Review of Ecology, Evolution, and
#' Systematics, 45, 297-324.
#'
#' @seealso [alphaDiversity()], [gammaDiversity()], [spaccBeta()] for
#'   spatial beta diversity accumulation
#'
#' @examples
#' # Simulated community data
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' # Partition diversity
#' part <- diversityPartition(species, q = c(0, 1, 2))
#' print(part)
#'
#' # Beta near 1 = homogeneous region
#' # Beta near n_sites = heterogeneous region
#'
#' @export
diversityPartition <- function(x, q = c(0, 1, 2), weights = "equal", coords = NULL) {
  x <- as.matrix(x)
  n_sites <- nrow(x)
  n_species <- ncol(x)

  if (!is.null(coords)) {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    stopifnot("x and coords must have same number of rows" = nrow(x) == nrow(coords))
  }

  # Compute per-site alpha
  alpha_per_site <- alphaDiversity(x, q)
  # Handle weights
  if (is.character(weights) && weights == "equal") {
    w <- rep(1 / n_sites, n_sites)
  } else if (is.character(weights) && weights == "proportional") {
    site_totals <- rowSums(x)
    w <- site_totals / sum(site_totals)
  } else if (is.numeric(weights)) {
    stopifnot(length(weights) == n_sites)
    w <- weights / sum(weights)
  } else {
    stop("weights must be 'equal', 'proportional', or a numeric vector")
  }
  # Weighted mean alpha (using generalized mean for Hill numbers)
  alpha <- sapply(seq_along(q), function(j) {
    if (q[j] == 1) {
      # Geometric mean for q=1
      exp(sum(w * log(alpha_per_site[[j]] + 1e-10)))
    } else {
      # Power mean for other q
      (sum(w * alpha_per_site[[j]]^(1 - q[j])))^(1 / (1 - q[j]))
    }
  })
  names(alpha) <- paste0("q", q)
  # Gamma diversity
  gamma <- gammaDiversity(x, q)
  # Beta = gamma / alpha (multiplicative partitioning)
  beta <- gamma / alpha
  names(beta) <- paste0("q", q)
  structure(
    list(
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      alpha_per_site = alpha_per_site,
      q = q,
      weights = w,
      coords = coords,
      n_sites = n_sites,
      n_species = n_species,
      call = match.call()
    ),
    class = "spacc_partition"
  )
}
#' @export
print.spacc_partition <- function(x, ...) {
  cat("Alpha-Beta-Gamma Diversity Partitioning\n")
  cat(sprintf("%d sites, %d species\n\n", x$n_sites, x$n_species))
  df <- data.frame(
    q = x$q,
    alpha = round(x$alpha, 2),
    beta = round(x$beta, 2),
    gamma = round(x$gamma, 2)
  )
  print(df, row.names = FALSE)
  cat("\nInterpretation:\n")
  cat("  Alpha = mean effective species per site\n")
  cat("  Beta  = effective number of communities (1 to n_sites)\n")
  cat("  Gamma = regional effective species (gamma = alpha x beta)\n")
  invisible(x)
}
#' @export
summary.spacc_partition <- function(object, ...) {
  data.frame(
    q = object$q,
    alpha = object$alpha,
    beta = object$beta,
    gamma = object$gamma,
    n_sites = object$n_sites,
    beta_normalized = (object$beta - 1) / (object$n_sites - 1)
  )
}


# spacc_alpha S3 methods ---------------------------------------------------

#' @export
print.spacc_alpha <- function(x, ...) {
  cat(sprintf("spacc alpha diversity: %d sites, %d species\n",
              x$n_sites, x$n_species))
  cat(sprintf("Orders (q): %s\n", paste(x$q, collapse = ", ")))
  if (!is.null(x$coords)) cat("Coordinates: available\n")
  invisible(x)
}


#' @export
summary.spacc_alpha <- function(object, ...) {
  q_cols <- paste0("q", object$q)
  vals <- object$values[, q_cols, drop = FALSE]

  data.frame(
    q = object$q,
    mean = colMeans(vals),
    sd = apply(vals, 2, stats::sd),
    min = apply(vals, 2, min),
    max = apply(vals, 2, max)
  )
}


#' @export
plot.spacc_alpha <- function(x, type = c("map", "histogram"), q = NULL,
                              point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (is.null(q)) q <- x$q[1]
  q_col <- paste0("q", q)
  stopifnot("q value not found" = q_col %in% names(x$values))

  if (type == "map") {
    if (is.null(x$coords)) stop("No coordinates available. Rerun alphaDiversity() with coords.")
    plot_spatial_map(x$values, q_col,
                     title = sprintf("Alpha diversity (q = %s)", q),
                     subtitle = sprintf("%d sites", x$n_sites),
                     point_size = point_size, palette = palette)
  } else {
    check_suggests("ggplot2")
    ggplot2::ggplot(x$values, ggplot2::aes(x = .data[[q_col]])) +
      ggplot2::geom_histogram(bins = 30, fill = "#4CAF50", color = "white", alpha = 0.8) +
      ggplot2::labs(x = sprintf("Alpha diversity (q = %s)", q), y = "Count",
                    title = sprintf("Distribution: alpha diversity (q = %s)", q)) +
      ggplot2::theme_minimal(base_size = 12)
  }
}


#' @export
as_sf.spacc_alpha <- function(x, crs = NULL) {
  if (is.null(x$coords)) stop("No coordinates available. Rerun alphaDiversity() with coords.")
  as_sf_from_df(x$values, crs = crs)
}


# spacc_partition spatial methods ------------------------------------------

#' @export
plot.spacc_partition <- function(x, type = c("map", "bar"), q = NULL,
                                  point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$coords)) stop("No coordinates available. Rerun diversityPartition() with coords.")
    if (is.null(q)) q <- x$q[1]
    q_col <- paste0("q", q)
    stopifnot("q value not found" = q_col %in% names(x$alpha_per_site))

    df <- x$alpha_per_site
    df$x <- x$coords$x
    df$y <- x$coords$y

    plot_spatial_map(df, q_col,
                     title = sprintf("Per-site alpha diversity (q = %s)", q),
                     subtitle = sprintf("%d sites", x$n_sites),
                     point_size = point_size, palette = palette)
  } else {
    check_suggests("ggplot2")
    summ <- summary(x)
    df <- data.frame(
      q = rep(paste0("q=", summ$q), 3),
      component = rep(c("Alpha", "Beta", "Gamma"), each = nrow(summ)),
      value = c(summ$alpha, summ$beta, summ$gamma)
    )
    df$component <- factor(df$component, levels = c("Gamma", "Alpha", "Beta"))

    ggplot2::ggplot(df, ggplot2::aes(x = q, y = value, fill = component)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Diversity order", y = "Effective species", fill = "Component",
                    title = "Diversity Partitioning") +
      ggplot2::theme_minimal(base_size = 12)
  }
}


#' @export
as_sf.spacc_partition <- function(x, crs = NULL) {
  if (is.null(x$coords)) stop("No coordinates available. Rerun diversityPartition() with coords.")
  df <- x$alpha_per_site
  df$x <- x$coords$x
  df$y <- x$coords$y
  as_sf_from_df(df, crs = crs)
}


# ============================================================================
# DIVERSITY PROFILES
# ============================================================================

#' Diversity Profile
#'
#' Compute Hill numbers across a continuous range of diversity orders (q),
#' producing a diversity profile for each site and the regional pool.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity to evaluate. Default
#'   `seq(0, 3, by = 0.1)`.
#' @param type Character. What to compute: `"per_site"` (per-site profiles),
#'   `"regional"` (pooled gamma), or `"both"` (default).
#' @param coords Optional data.frame with columns `x` and `y` for spatial
#'   mapping. When provided, enables `plot(type = "map")`.
#'
#' @return An object of class `spacc_profile` containing:
#'   \item{per_site}{Matrix of per-site diversity (sites x q values), or `NULL`}
#'   \item{regional}{Named numeric vector of gamma diversity per q, or `NULL`}
#'   \item{q}{Vector of diversity orders used}
#'   \item{coords}{Coordinates if provided}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Number of species}
#'
#' @details
#' A diversity profile plots effective number of species as a function of the
#' order q. The key property is that Hill numbers are non-increasing in q:
#' \eqn{D_q \ge D_{q'}} for \eqn{q \le q'}.
#'
#' - q = 0: Species richness (insensitive to abundance)
#' - q = 1: Exponential of Shannon entropy (all species weighted equally by
#'   their proportional abundance)
#' - q = 2: Inverse Simpson concentration (emphasizes dominant species)
#' - q > 2: Increasingly dominated by common species
#'
#' @references
#' Leinster, T. & Cobbold, C.A. (2012). Measuring diversity: the importance
#' of species similarity. Ecology, 93, 477-489.
#'
#' Chao, A., Chiu, C.H. & Jost, L. (2014). Unifying species diversity,
#' phylogenetic diversity, functional diversity, and related similarity and
#' differentiation measures through Hill numbers. Annual Review of Ecology,
#' Evolution, and Systematics, 45, 297-324.
#'
#' @seealso [alphaDiversity()] for per-site values at specific q,
#'   [gammaDiversity()] for regional diversity, [evenness()] for evenness
#'   profiles
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' prof <- diversityProfile(species)
#' print(prof)
#'
#' \donttest{
#' plot(prof)
#' }
#'
#' @export
diversityProfile <- function(x, q = seq(0, 3, by = 0.1),
                              type = c("both", "per_site", "regional"),
                              coords = NULL) {
  type <- match.arg(type)
  x <- as.matrix(x)
  n_sites <- nrow(x)
  n_species <- ncol(x)

  if (!is.null(coords)) {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    stopifnot("x and coords must have same number of rows" = nrow(x) == nrow(coords))
  }

  per_site <- NULL
  regional <- NULL

  if (type %in% c("both", "per_site")) {
    per_site <- matrix(NA, nrow = n_sites, ncol = length(q))
    colnames(per_site) <- paste0("q", q)
    for (i in seq_len(n_sites)) {
      abundances <- as.numeric(x[i, ])
      for (j in seq_along(q)) {
        per_site[i, j] <- calc_hill_number(abundances, q[j])
      }
    }
  }

  if (type %in% c("both", "regional")) {
    pooled <- colSums(x)
    regional <- sapply(q, function(qi) calc_hill_number(pooled, qi))
    names(regional) <- paste0("q", q)
  }

  structure(
    list(
      per_site = per_site,
      regional = regional,
      q = q,
      type = type,
      coords = coords,
      n_sites = n_sites,
      n_species = n_species,
      call = match.call()
    ),
    class = "spacc_profile"
  )
}


#' @export
print.spacc_profile <- function(x, ...) {
  cat(sprintf("Diversity profile: %d sites, %d species\n", x$n_sites, x$n_species))
  cat(sprintf("q range: [%.1f, %.1f] (%d values)\n", min(x$q), max(x$q), length(x$q)))
  if (!is.null(x$per_site)) {
    cat(sprintf("Per-site: mean D_0 = %.1f, mean D_1 = %.1f, mean D_2 = %.1f\n",
                mean(x$per_site[, 1]),
                mean(x$per_site[, which.min(abs(x$q - 1))]),
                mean(x$per_site[, which.min(abs(x$q - 2))])))
  }
  if (!is.null(x$regional)) {
    cat(sprintf("Regional: D_0 = %.1f, D_1 = %.1f, D_2 = %.1f\n",
                x$regional[1],
                x$regional[which.min(abs(x$q - 1))],
                x$regional[which.min(abs(x$q - 2))]))
  }
  invisible(x)
}


#' @export
summary.spacc_profile <- function(object, ...) {
  df <- data.frame(q = object$q)
  if (!is.null(object$per_site)) {
    df$mean_alpha <- colMeans(object$per_site)
    df$sd_alpha <- apply(object$per_site, 2, stats::sd)
    df$min_alpha <- apply(object$per_site, 2, min)
    df$max_alpha <- apply(object$per_site, 2, max)
  }
  if (!is.null(object$regional)) {
    df$gamma <- as.numeric(object$regional)
  }
  df
}


#' @export
as.data.frame.spacc_profile <- function(x, row.names = NULL, optional = FALSE, ...) {
  summary(x)
}


#' @export
as_sf.spacc_profile <- function(x, crs = NULL) {
  if (is.null(x$coords)) stop("No coordinates available. Rerun diversityProfile() with coords.")
  if (is.null(x$per_site)) stop("No per-site data. Rerun diversityProfile() with type = 'both' or 'per_site'.")
  df <- as.data.frame(x$per_site)
  df$x <- x$coords$x
  df$y <- x$coords$y
  as_sf_from_df(df, crs = crs)
}


#' @export
plot.spacc_profile <- function(x, type = c("profile", "map"),
                                show_sites = FALSE, ci = TRUE,
                                q_map = 0,
                                point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$coords)) stop("No coordinates available. Rerun diversityProfile() with coords.")
    if (is.null(x$per_site)) stop("No per-site data. Rerun diversityProfile() with type = 'both' or 'per_site'.")

    q_idx <- which.min(abs(x$q - q_map))
    q_col <- paste0("q", x$q[q_idx])
    df <- data.frame(x = x$coords$x, y = x$coords$y)
    df[[q_col]] <- x$per_site[, q_idx]

    return(plot_spatial_map(df, q_col,
                            title = sprintf("Diversity (q = %s)", x$q[q_idx]),
                            subtitle = sprintf("%d sites", x$n_sites),
                            point_size = point_size, palette = palette))
  }

  check_suggests("ggplot2")

  summ <- summary(x)
  p <- ggplot2::ggplot()

  if (!is.null(x$per_site)) {
    if (show_sites) {
      site_df <- data.frame(
        q = rep(x$q, each = x$n_sites),
        diversity = as.vector(x$per_site),
        site = rep(seq_len(x$n_sites), length(x$q))
      )
      p <- p + ggplot2::geom_line(
        data = site_df,
        ggplot2::aes(x = q, y = diversity, group = site),
        alpha = 0.1, color = "grey50"
      )
    }

    if (ci && "sd_alpha" %in% names(summ)) {
      p <- p + ggplot2::geom_ribbon(
        data = summ,
        ggplot2::aes(x = q, ymin = mean_alpha - sd_alpha,
                     ymax = mean_alpha + sd_alpha),
        alpha = 0.3, fill = "#4CAF50"
      )
    }

    p <- p + ggplot2::geom_line(
      data = summ,
      ggplot2::aes(x = q, y = mean_alpha),
      linewidth = 1, color = "#2E7D32"
    )
  }

  if (!is.null(x$regional)) {
    p <- p + ggplot2::geom_line(
      data = summ,
      ggplot2::aes(x = q, y = gamma),
      linewidth = 1, color = "#1565C0", linetype = "dashed"
    )
  }

  p +
    ggplot2::labs(
      x = "Diversity order (q)",
      y = "Effective number of species",
      title = "Diversity Profile"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
