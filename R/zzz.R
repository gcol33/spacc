#' @keywords internal
"_PACKAGE"

#' @useDynLib spacc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats quantile sd coef confint predict nls nls.control AIC BIC anova dist median runif lm approx
#' @importFrom parallel detectCores
NULL

# Global variable declarations for ggplot2 NSE and C++ functions
utils::globalVariables(c(
  # ggplot2 aesthetics
  "sites", "species", "seed", "lower", "upper", "y", "mean", "area",
  "value", "component", "beta_total", "beta_total_lower", "beta_total_upper",
  "mean_richness", "richness_lower", "richness_upper",
  "distance", "i", "j", "metric", "q_label", "group", "n", "expected",
  ".data",
  # C++ functions from RcppExports.R
  "calc_chiu_coverage", "calc_coverage", "calc_evo_heritage", "calc_faith_pd",
  "calc_hill_at_coverage", "calc_hill_at_size", "calc_k_facets", "calc_fdis", "calc_fric_approx",
  "calc_hill_number", "calc_mntd", "calc_mpd", "calc_rao", "calc_rao_traits",
  "cpp_beta_order_parallel",
  "cpp_cone_order", "cpp_distance_matrix",
  "cpp_focal_order", "cpp_func_order_parallel", "cpp_gaussian_order",
  "cpp_kncn_order", "cpp_kncn_tree_order",
  "cpp_nn_walk_order", "cpp_nn_walk_tree_order",
  "cpp_order_coverage_parallel", "cpp_order_hill_beta_parallel",
  "cpp_order_hill_coverage_parallel", "cpp_order_hill_parallel",
  "cpp_order_parallel", "cpp_phylo_order_parallel", "cpp_radius_parallel",
  "interpolate_at_coverage",
  # ggplot2 aesthetics for new features
  "mean_area", "mean_diversity", "mean_endemism", "mean_richness",
  "endemism_lower", "endemism_upper", "type", "predicted",
  "coverage", "richness", "target_coverage", "q_label",
  # Feature 1: Richness estimators
  "estimator", "estimate", "label",
  # Feature 2-3: Diversity and evenness profiles
  "diversity", "regional", "mean_alpha", "sd_alpha", "gamma",
  "evenness", "E_q", "pielou", "simpson_even", "site",
  # Feature 4: Beta distance-decay
  "dissimilarity", "site_i", "site_j", "predicted_decay", "bin_mid", "bin",
  # Feature 5: Null models / SES
  "ses", "null_mean", "null_sd", "p_value", "position", "observed",
  "significance",
  # Hill + coverage, Hill beta, MEM
  "mean_hill", "hill_lower", "hill_upper", "mean_coverage",
  "coverage_lower", "coverage_upper",
  "mean_gamma", "mean_alpha", "mean_beta",
  "gamma_lower", "gamma_upper", "alpha_lower", "alpha_upper",
  "beta_lower", "beta_upper",
  "eigenvalue", "moran_i", "vector", "r_squared", "score",
  # compareModels plot aesthetics
  "converged",
  "model",
  # extrapolate / extrapolateArea prediction bands
  "fit", "n_sub", "richness_sd"
))


# Null-coalescing helper (base R gained %||% in 4.4.0; define for R >= 3.5)
`%||%` <- function(a, b) if (is.null(a)) b else a


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
    message("i ", msg)
  }
}

cli_success <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_success(msg)
  } else {
    message("v ", msg)
  }
}

# Internal ggplot2 theme: minimal with no grid, transparent background
spacc_theme <- function(base_size = 19) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}
