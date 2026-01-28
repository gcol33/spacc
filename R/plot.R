#' Plot Spatial SAC
#'
#' Create a ggplot2 visualization of species accumulation curves.
#' For grouped `spacc` objects, curves are overlaid with different colors.
#'
#' @param x A `spacc` or `summary.spacc` object.
#' @param ci Logical. Show confidence interval ribbon? Default `TRUE`.
#' @param ci_level Numeric. Confidence level for interval. Default 0.95.
#' @param ci_alpha Numeric. Transparency of CI ribbon. Default 0.3.
#' @param show_seeds Logical. Show individual seed curves? Default `FALSE`.
#'   Only available for ungrouped objects.
#' @param saturation Logical. Mark saturation point? Default `FALSE`.
#'   Only available for ungrouped objects.
#' @param saturation_level Numeric. Proportion for saturation. Default 0.9.
#' @param facet Logical. Use faceted panels for grouped objects? Default `FALSE`
#'   (overlay with colors). Ignored for ungrouped objects.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' sac <- spacc(species, coords)
#' plot(sac)
#'
#' # Grouped
#' sac_grouped <- spacc(species, coords, groups = status)
#' plot(sac_grouped)
#' plot(sac_grouped, facet = TRUE)
#' }
#'
#' @export
plot.spacc <- function(x,
                       ci = TRUE,
                       ci_level = 0.95,
                       ci_alpha = 0.3,
                       show_seeds = FALSE,
                       saturation = FALSE,
                       saturation_level = 0.9,
                       facet = FALSE,
                       ...) {

  check_suggests("ggplot2")

  # Grouped spacc: overlay or facet

  if (is_grouped(x)) {
    summaries <- summary(x, ci_level = ci_level)

    df_list <- lapply(names(summaries), function(g) {
      s <- summaries[[g]]
      # Skip empty groups (0 sites)
      if (length(s$sites) == 0) return(NULL)
      data.frame(
        sites = s$sites,
        mean = s$mean,
        lower = s$lower,
        upper = s$upper,
        group = g
      )
    })
    # Remove NULL entries from empty groups
    df_list <- df_list[!vapply(df_list, is.null, logical(1))]

    df <- do.call(rbind, df_list)
    df$group <- factor(df$group, levels = x$group_names)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = sites, y = mean, color = group, fill = group))

    if (ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = ci_alpha, color = NA
      )
    }

    p <- p +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(
        x = "Sites sampled",
        y = "Cumulative species",
        title = "Species Accumulation Curves",
        color = "Group",
        fill = "Group"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    if (facet) {
      p <- p + ggplot2::facet_wrap(~ group, scales = "free_x")
    }

    return(p)
  }

  # Ungrouped spacc
  summ <- summary(x, saturation_threshold = saturation_level, ci_level = ci_level)

  df <- data.frame(
    sites = summ$sites,
    mean = summ$mean,
    lower = summ$lower,
    upper = summ$upper
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sites, y = mean))

  # Individual curves
 if (show_seeds) {
    curves_long <- data.frame(
      sites = rep(seq_len(ncol(x$curves)), each = nrow(x$curves)),
      species = as.vector(t(x$curves)),
      seed = rep(seq_len(nrow(x$curves)), ncol(x$curves))
    )
    p <- p + ggplot2::geom_line(
      data = curves_long,
      ggplot2::aes(x = sites, y = species, group = seed),
      alpha = 0.1, color = "grey50"
    )
  }

  # CI ribbon
  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }

  # Mean line
  p <- p + ggplot2::geom_line(linewidth = 1, color = "#2E7D32")

  # Saturation point
  if (saturation) {
    sat_df <- data.frame(
      x = summ$saturation_point,
      y = summ$mean[summ$saturation_point]
    )
    p <- p +
      ggplot2::geom_vline(
        xintercept = summ$saturation_point,
        linetype = "dashed", color = "#FF9800", linewidth = 0.8
      ) +
      ggplot2::geom_point(
        data = sat_df,
        ggplot2::aes(x = x, y = y),
        size = 3, color = "#FF9800"
      )
  }

  p <- p +
    ggplot2::labs(
      x = "Sites sampled",
      y = "Cumulative species",
      title = "Species Accumulation Curve",
      subtitle = sprintf("%s, %d seeds, %.0f%% CI",
                        x$method, x$n_seeds, ci_level * 100)
    ) +
    ggplot2::theme_minimal(base_size = 12)

  p
}
