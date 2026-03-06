#' Evenness Profiles
#'
#' Compute species evenness across sites using Hill-based, Pielou, or
#' Simpson evenness measures.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity for Hill evenness.
#'   Default `seq(0.1, 3, by = 0.1)`. Note: q = 0 is excluded by default
#'   because Hill evenness at q = 0 is trivially S/S = 1.
#' @param type Character. Evenness type: `"hill"` (Hill evenness E_q = D_q / D_0,
#'   default), `"pielou"` (Pielou's J = log(D_1) / log(S)), or
#'   `"simpson"` (Simpson evenness = (1/D_2) / S).
#' @param coords Optional data.frame with columns `x` and `y` for spatial
#'   mapping. When provided, enables `plot(type = "map")`.
#'
#' @return An object of class `spacc_evenness` containing:
#'   \item{per_site}{Matrix or vector of per-site evenness values}
#'   \item{regional}{Regional (pooled) evenness}
#'   \item{q}{Orders used (for Hill type)}
#'   \item{type}{Evenness type}
#'   \item{coords}{Coordinates if provided}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Number of species}
#'
#' @details
#' All evenness measures are bounded in \eqn{[0, 1]}:
#' - 0 = maximally uneven (one dominant species)
#' - 1 = perfectly even (all species equally abundant)
#'
#' **Hill evenness** (Jost 2010):
#' \deqn{E_q = D_q / D_0}
#' This divides the effective number of species at order q by species richness.
#'
#' **Pielou's J** (Pielou 1966):
#' \deqn{J = \frac{\log(D_1)}{\log(S)} = \frac{H'}{\log(S)}}
#'
#' **Simpson evenness**:
#' \deqn{E_{1/D} = \frac{1}{D_2 \cdot S}}
#'
#' @references
#' Jost, L. (2010). The relation between evenness and diversity. Diversity, 2,
#' 207-232.
#'
#' Pielou, E.C. (1966). The measurement of diversity in different types of
#' biological collections. Journal of Theoretical Biology, 13, 131-144.
#'
#' @seealso [diversityProfile()] for Hill number profiles,
#'   [alphaDiversity()] for raw diversity values
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' # Hill evenness profile
#' ev <- evenness(species)
#' print(ev)
#'
#' # Pielou's J
#' ev_j <- evenness(species, type = "pielou")
#' print(ev_j)
#'
#' @export
evenness <- function(x, q = seq(0.1, 3, by = 0.1),
                     type = c("hill", "pielou", "simpson"),
                     coords = NULL) {
  type <- match.arg(type)
  x <- as.matrix(x)
  n_sites <- nrow(x)
  n_species <- ncol(x)

  if (!is.null(coords)) {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    stopifnot("x and coords must have same number of rows" = nrow(x) == nrow(coords))
  }

  # Compute richness (D_0) per site
  richness <- rowSums(x > 0)

  if (type == "hill") {
    # Hill evenness: E_q = D_q / D_0
    per_site <- matrix(NA, nrow = n_sites, ncol = length(q))
    colnames(per_site) <- paste0("q", q)

    for (i in seq_len(n_sites)) {
      abundances <- as.numeric(x[i, ])
      S <- richness[i]
      if (S <= 1) {
        per_site[i, ] <- 1  # Perfectly even by convention for S <= 1
      } else {
        for (j in seq_along(q)) {
          Dq <- calc_hill_number(abundances, q[j])
          per_site[i, j] <- Dq / S
        }
      }
    }

    # Regional
    pooled <- colSums(x)
    S_regional <- sum(pooled > 0)
    if (S_regional <= 1) {
      regional <- rep(1, length(q))
    } else {
      regional <- sapply(q, function(qi) calc_hill_number(pooled, qi)) / S_regional
    }
    names(regional) <- paste0("q", q)

  } else if (type == "pielou") {
    # Pielou's J = log(D_1) / log(S)
    per_site <- numeric(n_sites)
    for (i in seq_len(n_sites)) {
      abundances <- as.numeric(x[i, ])
      S <- richness[i]
      if (S <= 1) {
        per_site[i] <- 1
      } else {
        D1 <- calc_hill_number(abundances, 1)
        per_site[i] <- log(D1) / log(S)
      }
    }

    pooled <- colSums(x)
    S_regional <- sum(pooled > 0)
    if (S_regional <= 1) {
      regional <- 1
    } else {
      D1_reg <- calc_hill_number(pooled, 1)
      regional <- log(D1_reg) / log(S_regional)
    }
    q <- 1  # Pielou is effectively q=1

  } else {
    # Simpson evenness = (1/D_2) / S = 1 / (D_2 * S)
    per_site <- numeric(n_sites)
    for (i in seq_len(n_sites)) {
      abundances <- as.numeric(x[i, ])
      S <- richness[i]
      if (S <= 1) {
        per_site[i] <- 1
      } else {
        D2 <- calc_hill_number(abundances, 2)
        per_site[i] <- 1 / (D2 * S)
      }
    }

    pooled <- colSums(x)
    S_regional <- sum(pooled > 0)
    if (S_regional <= 1) {
      regional <- 1
    } else {
      D2_reg <- calc_hill_number(pooled, 2)
      regional <- 1 / (D2_reg * S_regional)
    }
    q <- 2  # Simpson is effectively q=2
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
    class = "spacc_evenness"
  )
}


#' @export
print.spacc_evenness <- function(x, ...) {
  cat(sprintf("Evenness (%s): %d sites, %d species\n", x$type, x$n_sites, x$n_species))
  if (x$type == "hill") {
    cat(sprintf("q range: [%.1f, %.1f]\n", min(x$q), max(x$q)))
    # Report at a few key q values
    q1_idx <- which.min(abs(x$q - 1))
    q2_idx <- which.min(abs(x$q - 2))
    cat(sprintf("Mean evenness: E_1 = %.3f, E_2 = %.3f\n",
                mean(x$per_site[, q1_idx]), mean(x$per_site[, q2_idx])))
    cat(sprintf("Regional evenness: E_1 = %.3f, E_2 = %.3f\n",
                x$regional[q1_idx], x$regional[q2_idx]))
  } else if (x$type == "pielou") {
    cat(sprintf("Mean Pielou's J: %.3f\n", mean(x$per_site)))
    cat(sprintf("Regional J: %.3f\n", x$regional))
  } else {
    cat(sprintf("Mean Simpson evenness: %.3f\n", mean(x$per_site)))
    cat(sprintf("Regional Simpson evenness: %.3f\n", x$regional))
  }
  invisible(x)
}


#' @export
summary.spacc_evenness <- function(object, ...) {
  if (object$type == "hill") {
    data.frame(
      q = object$q,
      mean = colMeans(object$per_site),
      sd = apply(object$per_site, 2, stats::sd),
      min = apply(object$per_site, 2, min),
      max = apply(object$per_site, 2, max),
      regional = as.numeric(object$regional)
    )
  } else {
    data.frame(
      type = object$type,
      mean = mean(object$per_site),
      sd = stats::sd(object$per_site),
      min = min(object$per_site),
      max = max(object$per_site),
      regional = object$regional
    )
  }
}


#' @export
as.data.frame.spacc_evenness <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (x$type == "hill") {
    data.frame(
      q = x$q,
      mean = colMeans(x$per_site),
      sd = apply(x$per_site, 2, stats::sd),
      regional = as.numeric(x$regional)
    )
  } else {
    data.frame(
      site = seq_len(x$n_sites),
      evenness = x$per_site
    )
  }
}


#' @export
as_sf.spacc_evenness <- function(x, crs = NULL) {
  if (is.null(x$coords)) stop("No coordinates available. Rerun evenness() with coords.")

  if (x$type == "hill") {
    df <- as.data.frame(x$per_site)
    colnames(df) <- paste0("E_", colnames(df))
  } else {
    df <- data.frame(evenness = x$per_site)
    names(df) <- paste0(x$type, "_evenness")
  }
  df$x <- x$coords$x
  df$y <- x$coords$y
  as_sf_from_df(df, crs = crs)
}


#' @export
plot.spacc_evenness <- function(x, type = c("profile", "map", "histogram"),
                                 q_map = 1,
                                 point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$coords)) stop("No coordinates available. Rerun evenness() with coords.")

    if (x$type == "hill") {
      q_idx <- which.min(abs(x$q - q_map))
      q_col <- paste0("E_q", x$q[q_idx])
      df <- data.frame(x = x$coords$x, y = x$coords$y)
      df[[q_col]] <- x$per_site[, q_idx]
      return(plot_spatial_map(df, q_col,
                              title = sprintf("Hill evenness (q = %s)", x$q[q_idx]),
                              subtitle = sprintf("%d sites", x$n_sites),
                              point_size = point_size, palette = palette))
    } else {
      q_col <- paste0(x$type, "_evenness")
      df <- data.frame(x = x$coords$x, y = x$coords$y)
      df[[q_col]] <- x$per_site
      return(plot_spatial_map(df, q_col,
                              title = sprintf("%s evenness", tools::toTitleCase(x$type)),
                              subtitle = sprintf("%d sites", x$n_sites),
                              point_size = point_size, palette = palette))
    }
  }

  check_suggests("ggplot2")

  if (type == "histogram") {
    if (x$type == "hill") {
      q_idx <- which.min(abs(x$q - q_map))
      vals <- x$per_site[, q_idx]
      label <- sprintf("Hill evenness (q = %s)", x$q[q_idx])
    } else {
      vals <- x$per_site
      label <- sprintf("%s evenness", tools::toTitleCase(x$type))
    }

    df <- data.frame(evenness = vals)
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data[["evenness"]])) +
        ggplot2::geom_histogram(bins = 30, fill = "#4CAF50", color = "white", alpha = 0.8) +
        ggplot2::labs(x = label, y = "Count", title = paste("Distribution:", label)) +
        spacc_theme()
    )
  }

  # Profile plot (only meaningful for Hill evenness)
  if (x$type != "hill") {
    stop("Profile plot only available for Hill evenness (type = 'hill'). Use type = 'histogram' instead.")
  }

  summ <- summary(x)

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = q, y = mean)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean - sd, ymax = pmin(mean + sd, 1)),
      alpha = 0.3, fill = "#4CAF50"
    ) +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::geom_line(
      ggplot2::aes(y = regional),
      linewidth = 1, color = "#1565C0", linetype = "dashed"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Diversity order (q)",
      y = "Evenness (E_q = D_q / D_0)",
      title = "Evenness Profile"
    ) +
    spacc_theme()

  p
}
