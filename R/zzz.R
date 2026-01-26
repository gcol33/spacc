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
  "cpp_beta_knn_parallel", "cpp_func_knn_parallel", "cpp_kncn_metrics_parallel",
 "cpp_knn_coverage_parallel", "cpp_knn_hill_parallel", "cpp_knn_metrics_parallel",
  "cpp_phylo_knn_parallel", "interpolate_at_coverage"
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
