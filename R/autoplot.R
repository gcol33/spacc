#' Autoplot Methods for spacc Objects
#'
#' `ggplot2::autoplot()` methods for spacc objects. Each one calls the
#' corresponding `plot()` method, which already returns a ggplot2 object.
#'
#' @name autoplot.spacc
#' @param object A spacc object.
#' @param ... Additional arguments passed to the corresponding `plot()` method.
#' @return A ggplot2 object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   coords <- data.frame(x = runif(30), y = runif(30))
#'   species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)
#'   sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
#'   ggplot2::autoplot(sac)
#' }
#' }
#'
#' @rawNamespace S3method(ggplot2::autoplot, spacc)
#' S3method(ggplot2::autoplot, spacc_beta)
#' S3method(ggplot2::autoplot, spacc_completeness)
#' S3method(ggplot2::autoplot, spacc_coverage)
#' S3method(ggplot2::autoplot, spacc_decay)
#' S3method(ggplot2::autoplot, spacc_endemism)
#' S3method(ggplot2::autoplot, spacc_fit)
#' S3method(ggplot2::autoplot, spacc_hill)
#' S3method(ggplot2::autoplot, spacc_hill_beta)
#' S3method(ggplot2::autoplot, spacc_hill_coverage)
#' S3method(ggplot2::autoplot, spacc_mem)
#' S3method(ggplot2::autoplot, spacc_mem_partition)
#' S3method(ggplot2::autoplot, spacc_metrics)
#' S3method(ggplot2::autoplot, spacc_model_compare)
#' S3method(ggplot2::autoplot, spacc_rare)
#' S3method(ggplot2::autoplot, spacc_redundancy)
#' S3method(ggplot2::autoplot, spacc_ses)
#' S3method(ggplot2::autoplot, spacc_zeta)
NULL

.autoplot_via_plot <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}

# One autoplot method per class, all delegating to plot(). Adding a class
# means adding its name here and an S3method line in the block above.
.autoplot_classes <- c(
  "spacc", "spacc_beta", "spacc_completeness", "spacc_coverage", "spacc_decay",
  "spacc_endemism", "spacc_fit", "spacc_hill", "spacc_hill_beta",
  "spacc_hill_coverage", "spacc_mem", "spacc_mem_partition", "spacc_metrics",
  "spacc_model_compare", "spacc_rare", "spacc_redundancy", "spacc_ses", "spacc_zeta"
)
for (.cls in .autoplot_classes) {
  assign(paste0("autoplot.", .cls), .autoplot_via_plot)
}
rm(.cls)
