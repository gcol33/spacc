#' Summarize Spatial SAC Results
#'
#' Compute summary statistics from multiple accumulation curves.
#'
#' @param object A `spatial_sac` object from [fast_spatial_sac()]
#' @param saturation_threshold Numeric. Proportion of maximum species to define
#'   saturation point. Default 0.9 (90%).
#' @param ci_level Numeric. Confidence interval level. Default 0.95.
#' @param ... Additional arguments (ignored)
#'
#' @return A `sac_summary` object containing:
#'   \item{mean_curve}{Mean accumulation across seeds}
#'   \item{lower}{Lower CI bound}
#'   \item{upper}{Upper CI bound}
#'   \item{saturation_point}{Mean number of sites to reach threshold}
#'   \item{cv_midpoint}{Coefficient of variation at midpoint}
#'
#' @export
sac_summary <- function(object, saturation_threshold = 0.9, ci_level = 0.95, ...) {
 UseMethod("sac_summary")
}

#' @export
sac_summary.spatial_sac <- function(object, saturation_threshold = 0.9, ci_level = 0.95, ...) {

 curves <- object$curves
 n_sites <- object$n_sites

 alpha <- (1 - ci_level) / 2

 # Compute statistics per site position
 mean_curve <- colMeans(curves)
 lower <- apply(curves, 2, quantile, probs = alpha)
 upper <- apply(curves, 2, quantile, probs = 1 - alpha)
 sd_curve <- apply(curves, 2, sd)

 # Saturation point: first site where mean reaches threshold * max
 max_species <- mean_curve[n_sites]
 threshold_value <- saturation_threshold * max_species
 saturation_point <- which(mean_curve >= threshold_value)[1]

 # CV at midpoint
 midpoint <- ceiling(n_sites / 2)
 cv_midpoint <- sd_curve[midpoint] / mean_curve[midpoint] * 100

 structure(
   list(
     sites = seq_len(n_sites),
     mean = mean_curve,
     lower = lower,
     upper = upper,
     sd = sd_curve,
     saturation_point = saturation_point,
     saturation_threshold = saturation_threshold,
     cv_midpoint = cv_midpoint,
     ci_level = ci_level,
     n_seeds = object$n_seeds,
     method = object$method
   ),
   class = "sac_summary"
 )
}

#' @export
print.sac_summary <- function(x, ...) {
 cat("Spatial SAC Summary\n")
 cat("-------------------\n")
 cat("Method:", x$method, "\n")
 cat("Seeds:", x$n_seeds, "\n")
 cat("Sites:", length(x$sites), "\n")
 cat("Final species:", round(x$mean[length(x$mean)], 1),
     sprintf("(%.0f%% CI: %.1f - %.1f)\n",
             x$ci_level * 100, x$lower[length(x$lower)], x$upper[length(x$upper)]))
 cat("Saturation (", x$saturation_threshold * 100, "%):", x$saturation_point, "sites\n", sep = "")
 cat("CV at midpoint:", round(x$cv_midpoint, 1), "%\n")
 invisible(x)
}
