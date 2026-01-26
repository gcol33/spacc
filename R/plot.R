#' Plot Spatial SAC
#'
#' Create a ggplot2 visualization of species accumulation curves.
#'
#' @param x A `spatial_sac` or `sac_summary` object
#' @param show_ci Logical. Show confidence interval ribbon? Default TRUE.
#' @param show_curves Logical. Show individual seed curves? Default FALSE.
#' @param ci_alpha Numeric. Transparency of CI ribbon. Default 0.3.
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_minimal
#' @export
plot.spatial_sac <- function(x, show_ci = TRUE, show_curves = FALSE, ci_alpha = 0.3, ...) {

 summ <- sac_summary(x)
 plot.sac_summary(summ, show_ci = show_ci, show_curves = show_curves,
                  ci_alpha = ci_alpha, curves = x$curves, ...)
}

#' @export
plot.sac_summary <- function(x, show_ci = TRUE, show_curves = FALSE,
                            ci_alpha = 0.3, curves = NULL, ...) {

 df <- data.frame(
   sites = x$sites,
   mean = x$mean,
   lower = x$lower,
   upper = x$upper
 )

 p <- ggplot2::ggplot(df, ggplot2::aes(x = sites, y = mean))

 # Optional: individual curves
 if (show_curves && !is.null(curves)) {
   curves_long <- data.frame(
     sites = rep(seq_len(ncol(curves)), each = nrow(curves)),
     species = as.vector(t(curves)),
     seed = rep(seq_len(nrow(curves)), ncol(curves))
   )
   p <- p + ggplot2::geom_line(
     data = curves_long,
     ggplot2::aes(x = sites, y = species, group = seed),
     alpha = 0.1, color = "grey50"
   )
 }

 # CI ribbon
 if (show_ci) {
   p <- p + ggplot2::geom_ribbon(
     ggplot2::aes(ymin = lower, ymax = upper),
     alpha = ci_alpha, fill = "#4CAF50"
   )
 }

 # Mean line
 p <- p +
   ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
   ggplot2::labs(
     x = "Sites sampled",
     y = "Cumulative species",
     title = "Spatial Species Accumulation Curve",
     subtitle = sprintf("%s method, %d seeds, %.0f%% CI",
                       x$method, x$n_seeds, x$ci_level * 100)
   ) +
   ggplot2::theme_minimal(base_size = 12)

 p
}
