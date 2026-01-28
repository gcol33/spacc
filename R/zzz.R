#' @keywords internal
"_PACKAGE"

#' @useDynLib spacc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats quantile sd coef confint predict nls nls.control AIC dist median runif lm
#' @importFrom parallel detectCores
NULL

# Global variable declarations for ggplot2 NSE and C++ functions
utils::globalVariables(c(
  # ggplot2 aesthetics
  "sites", "species", "seed", "lower", "upper", "y", "mean",
  "value", "component", "beta_total", "beta_total_lower", "beta_total_upper",
  "mean_richness", "richness_lower", "richness_upper",
  "distance", "i", "j", "metric", "q_label", "group", "n", "expected",
  ".data",
  # C++ functions from RcppExports.R
  "calc_coverage", "calc_faith_pd", "calc_fdis", "calc_fric_approx",
  "calc_hill_number", "calc_mntd", "calc_mpd",
  "cpp_beta_knn_parallel", "cpp_beta_knn_single",
  "cpp_collector_single", "cpp_cone_parallel", "cpp_cone_single",
  "cpp_distance_decay_parallel", "cpp_distance_decay_single",
  "cpp_distance_matrix",
  "cpp_func_knn_parallel", "cpp_func_knn_single",
  "cpp_gaussian_parallel", "cpp_gaussian_single",
  "cpp_kncn_kdtree_parallel", "cpp_kncn_kdtree_single",
  "cpp_kncn_metrics_parallel", "cpp_kncn_parallel", "cpp_kncn_single",
  "cpp_knn_coverage_parallel", "cpp_knn_coverage_single",
  "cpp_knn_hill_parallel", "cpp_knn_hill_single",
  "cpp_knn_kdtree_parallel", "cpp_knn_kdtree_parallel_seeds", "cpp_knn_kdtree_single",
  "cpp_knn_metrics_parallel", "cpp_knn_parallel", "cpp_knn_parallel_seeds", "cpp_knn_single",
  "cpp_phylo_knn_parallel", "cpp_phylo_knn_single",
  "cpp_radius_parallel", "cpp_radius_single",
  "cpp_random_parallel", "cpp_random_single",
  "cpp_wavefront_parallel", "cpp_wavefront_single",
  "interpolate_at_coverage",
  # ggplot2 aesthetics for new features
  "mean_area", "mean_diversity", "mean_endemism", "mean_richness",
  "endemism_lower", "endemism_upper", "type", "predicted",
  "coverage", "richness", "target_coverage", "q_label"
))


# Check for suggested packages
check_suggests <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' is required for this function. Install with: install.packages('%s')", pkg, pkg),
      call. = FALSE
    )
  }
}


# Safe cli messaging (falls back to cat if cli not available)
cli_info <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_info(msg)
  } else {
    cat("i", msg, "\n")
  }
}

cli_success <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_success(msg)
  } else {
    cat("v", msg, "\n")
  }
}
